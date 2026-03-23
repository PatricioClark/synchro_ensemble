!=================================================================
      PROGRAM MAIN
!=================================================================
! GHOST code: Geophysical High Order Suite for Turbulence
!
! Numerically integrates several fluid dynamics equations
! in 2 and 3 dimensions with periodic boundary conditions
! and external forcing. A pseudo-spectral method is used to
! compute spatial derivatives, while different time steppings
! can be usedd to evolve the system in the time domain.
!
! Notation: index 'i' is 'x'
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2003 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!
! References:
! Mininni PD, Rosenberg DL, Reddy R, Pouquet A.; P.Comp.37, 123 (2011)
! Rosenberg DL, Mininni PD, Reddy R, Pouquet A.: Atmosph.11, 178 (2020)
!=================================================================

!
! Definitions for conditional compilation
#include "ghost3D.h"

!
! Modules

      USE commtypes
      USE filefmt
      USE iovar
      USE fft
      USE threads
      USE offloading
      USE boxsize
      USE status
      USE gtimer
      USE ic_factory
      USE force_factory
      USE equation_factory
      USE stepper_factory
!     USE class_GPart

      IMPLICIT NONE

!
! Type encapsulating all state for a single nudged simulation member
      TYPE :: GNudgedSim
          TYPE(GStateComp),    ALLOCATABLE :: field(:),field_nxt(:)
          CLASS(EquationBase), ALLOCATABLE :: pde
          CLASS(GStepperBase), ALLOCATABLE :: stepper
          CLASS(icChain),      ALLOCATABLE :: iclist(:)
          INTEGER                          :: num_components
      END TYPE GNudgedSim

!
! Arrays for the field states, workspace, I/O, and PDE solver class

      TYPE(GStateComp),    ALLOCATABLE :: field(:),field_nxt(:),force(:),diff(:)
      TYPE(GWorkspace)                 :: workspace
      TYPE(IOPLAN)                     :: planio
      CLASS(EquationBase), ALLOCATABLE :: pde
      CLASS(GStepperBase), ALLOCATABLE :: stepper
      CLASS(icChain),      ALLOCATABLE :: iclist(:)
      CLASS(forceChain),   ALLOCATABLE :: forcemethod(:)

      TYPE(GNudgedSim), ALLOCATABLE :: ensemble(:)

!
! Auxiliary variables

      INTEGER :: idevice, iret, ncuda, ngcuda, ppn
      INTEGER :: num_components
      INTEGER :: num_realizations, kndg
      INTEGER :: ir, ip
      INTEGER :: t,timep,pstep,lgmult
      CHARACTER(LEN=8)   :: outlabel
      CHARACTER(LEN=128) :: odir_save
      INTEGER :: ihcpu1,ihcpu2
      INTEGER :: ihomp1,ihomp2
      INTEGER :: ihwtm1,ihwtm2
      LOGICAL :: bbenchexist
      NAMELIST / ensembleparams / num_realizations, kndg

!
! Initialization
! Initializes the MPI and I/O libraries
      CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

! Initializes the grid. This must be done early to have nx, ny, nz.
      CALL grid_init('parameter.inp')

! Initialization of offloading to GPUs using OpenMP (this is independent
! of CUDA initialization in systems with NVIDIA GPUs). GHOST
! assumes the number of MPI jobs in each node is equal to the
! number of GPUs available in the node. The user must ensure this
! condition is fulfilled.
#if defined(DO_HYBRIDoffl)
      CALL init_offload(myrank,numdev,hostdev,targetdev)
#endif

! Initialization of time integration parameters
      CALL status_init('parameter.inp')

! Read ensemble size from parameter file
      num_realizations = 1
      kndg      = 1
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form='formatted')
         READ(1,NML=ensembleparams)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(num_realizations,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kndg     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! I/O initialization
      CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,(/nx,ny,nz/),ksta,kend,planio)

! Now we can initialize the PDE method
     pde          = init_pdes_from_file(   'parameter.inp')
     CALL workspace%initialize_pool(NUMTMPREAL,NUMTMPCOMP)
     CALL pde%Solver_ctor('parameter.inp',workspace,planio)
     num_components = pde%state_size()
     CALL GState_alloc(field    , num_components)
     CALL GState_alloc(field_nxt, num_components)
     CALL GState_alloc(force    , num_components)
     CALL GState_alloc(diff     , num_components)
     ALLOCATE(icChain :: iclist(1))
     ALLOCATE(read_v  :: iclist(1)%ic)
     forcemethod  = init_forcing_from_file('parameter.inp',workspace)

     ALLOCATE(ensemble(num_realizations))
     DO ir = 1, SIZE(ensemble)
       ! BLOCK and MOVE_ALLOC are needed becauce the base class
       ! does not initilize the pointers to null ("=> null()" in declaration)
       ! This is a workaround
       BLOCK
         CLASS(EquationBase), ALLOCATABLE :: tmp_pde
         tmp_pde = init_pdes_from_file('parameter.inp')
         CALL MOVE_ALLOC(tmp_pde, ensemble(ir)%pde)
       END BLOCK
       CALL ensemble(ir)%pde%Solver_ctor('parameter.inp',workspace,planio)
       ensemble(ir)%num_components = ensemble(ir)%pde%state_size()
       CALL GState_alloc(ensemble(ir)%field    , ensemble(ir)%num_components)
       CALL GState_alloc(ensemble(ir)%field_nxt, ensemble(ir)%num_components)
       BLOCK
         CLASS(icChain),    ALLOCATABLE :: tmp_ic(:)
         tmp_ic  = init_ic_from_file(     'parameter.inp')
         CALL MOVE_ALLOC(tmp_ic,  ensemble(ir)%iclist)
       END BLOCK
     END DO

! Create ensemble output directories (odir_ensNNN, odir_difNNN)
     IF (myrank.eq.0) THEN
        DO ir = 1, SIZE(ensemble)
           WRITE(outlabel,'(A,I3.3)') '_ens', ir
           CALL EXECUTE_COMMAND_LINE('mkdir -p '//TRIM(odir)//TRIM(outlabel))
           WRITE(outlabel,'(A,I3.3)') '_dif', ir
           CALL EXECUTE_COMMAND_LINE('mkdir -p '//TRIM(odir)//TRIM(outlabel))
        END DO
     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

! Initialization of the numerical domain
     CALL box_init('parameter.inp')

! Initializes the FFT library. This must be done at
! this stage as it requires the variable "bench" to
! be properly initialized.
! Use FFTW_ESTIMATE or FFTW_MEASURE in short runs
! Use FFTW_PATIENT or FFTW_EXHAUSTIVE in long runs

      nth = 1
!$    nth = omp_get_max_threads()
#if !defined(DEF_GHOST_CUDA_)
!$    CALL fftp3d_init_threads(ierr)
#endif
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTInitHandle(ihcpu2,GT_CPUTIME)
         CALL GTInitHandle(ihomp2,GT_OMPTIME)
         CALL GTInitHandle(ihwtm2,GT_WTIME)
         CALL GTStart(ihcpu2)
         CALL GTStart(ihomp2)
         CALL GTStart(ihwtm2)
      ENDIF
      CALL fftp3d_create_plan(planrc,(/nx,ny,nz/),FFTW_REAL_TO_COMPLEX, &
                             FFTW_ESTIMATE)
      CALL fftp3d_create_plan(plancr,(/nx,ny,nz/),FFTW_COMPLEX_TO_REAL, &
                             FFTW_ESTIMATE)
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu2)
         CALL GTStop(ihomp2)
         CALL GTStop(ihwtm2)
      ENDIF

! Initial states
 IC : IF (stat.eq.0) THEN                 ! If stat=0 we start a new run
        ini  = 1
        sind = 0                          ! index for the spectrum
        tind = 0                          ! index for the binaries
!       pind = 0                          ! index for the particles
        timet = tstep
        timep = pstep
        timec = cstep
        times = sstep
        timef = 0
      ELSE                    ! If stat.ne.0 a previous run is continued
        ini = int((stat-1)*tstep) + 1
        tind = int(stat)
        sind = int(real(ini,kind=GP)/real(sstep,kind=GP)+1)
!       pind = int((stat-1)*lgmult+1)
        timet = 0
        timep = 0
        timec = int(modulo(float(ini-1),float(cstep)))
        times = int(modulo(float(ini-1),float(sstep)))
        timef = int(modulo(float(ini-1),float(fstep)))
      ENDIF IC
      CALL init_allstates(iclist,pde,field)
      CALL init_forcing(forcemethod,pde,force)

      DO ir = 1, SIZE(ensemble)
        CALL init_allstates(ensemble(ir)%iclist, ensemble(ir)%pde, ensemble(ir)%field)
        BLOCK
          CLASS(GStepperBase), ALLOCATABLE :: tmp_stp
          tmp_stp = build_stepper_from_file('parameter.inp',workspace,ensemble(ir)%pde)
          CALL MOVE_ALLOC(tmp_stp, ensemble(ir)%stepper)
        END BLOCK
      END DO

! Initial replace scales
     DO ir = 1, SIZE(ensemble)
        CALL replace_scales(field, ensemble(ir)%field, kndg)
     END DO

! Sets up the time stepper
      stepper = build_stepper_from_file('parameter.inp',workspace,pde)

! Time integration scheme starts here.
! If we are doing a benchmark, we measure
! cputime before starting.

      IF (bench.eq.1) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTInitHandle(ihcpu1,GT_CPUTIME)
         CALL GTInitHandle(ihomp1,GT_OMPTIME)
         CALL GTInitHandle(ihwtm1,GT_WTIME)
         ffttime  = 0.D00 ! re-inititialize fftp timers
         tratime  = 0.0D0
         comtime  = 0.D00
         tottime  = 0.0D0
#if defined(DEF_GHOST_CUDA_)
         memtime  = 0.0D0
         asstime  = 0.D00
#endif
         CALL GTStart(ihcpu1)
         CALL GTStart(ihomp1)
         CALL GTStart(ihwtm1)
      ENDIF

 RK : DO t = ini,step

! Every 'tstep' steps, stores the fields in binary files

         IF ((timet.eq.tstep).and.(bench.eq.0)) THEN
            timet = 0
            tind = tind+1
            CALL pde%write_states(field, planio)
            DO ir = 1, SIZE(ensemble)
               odir_save = odir
               WRITE(outlabel,'(A,I3.3)') '_ens', ir
               odir = TRIM(odir_save)//TRIM(outlabel)
               CALL ensemble(ir)%pde%write_states(ensemble(ir)%field, planio)
               DO ip = 1, num_components
                  diff(ip)%ccomp = ensemble(ir)%field(ip)%ccomp - field(ip)%ccomp
               END DO
               WRITE(outlabel,'(A,I3.3)') '_dif', ir
               odir = TRIM(odir_save)//TRIM(outlabel)
               CALL pde%write_states(diff, planio)
               odir = odir_save
            END DO
         ENDIF

! Every 'cstep' steps writes global quantities

         IF ((timec.eq.cstep).and.(bench.eq.0)) THEN
            timec = 0
            CALL hdcheck_ndg(field, force, t, dt, 0, 0, '_ref')
            DO ir = 1, SIZE(ensemble)
               WRITE(outlabel,'(A,I3.3,A)') '_ens', ir
               CALL hdcheck_ndg(ensemble(ir)%field, force, t, dt, 0, 0, outlabel)
               DO ip = 1, num_components
                  diff(ip)%ccomp = ensemble(ir)%field(ip)%ccomp - field(ip)%ccomp
               END DO
               WRITE(outlabel,'(A,I3.3,A)') '_dif', ir
               CALL hdcheck_ndg(diff, force, t, dt, 0, 0, outlabel)
            END DO
         ENDIF

! Every 'sstep' steps writes spectra

         IF ((times.eq.sstep).and.(bench.eq.0)) THEN
            times = 0
            sind = sind+1
            WRITE(ext, fmtext) sind
            CALL hdspectrum_ndg(field, '_ref')
            DO ir = 1, SIZE(ensemble)
               WRITE(outlabel,'(A,I3.3)') '_ens', ir
               CALL hdspectrum_ndg(ensemble(ir)%field, outlabel)
               DO ip = 1, num_components
                  diff(ip)%ccomp = ensemble(ir)%field(ip)%ccomp - field(ip)%ccomp
               END DO
               WRITE(outlabel,'(A,I3.3)') '_dif', ir
               CALL hdspectrum_ndg(diff, outlabel)
            END DO
         ENDIF

! Time evolution
         CALL update_forcing(forcemethod,pde,force)
         CALL stepper%step(time, field, force, dt, field_nxt)
    
         DO ir = 1, SIZE(ensemble)
            CALL replace_scales(field_nxt, ensemble(ir)%field, kndg)
            CALL ensemble(ir)%stepper%step(time, &
                                           ensemble(ir)%field, &
                                           force, &
                                           dt, &
                                           ensemble(ir)%field_nxt)
            ensemble(ir)%field = ensemble(ir)%field_nxt
         END DO

         field = field_nxt
         timet = timet+1
         timep = timep+1
         timec = timec+1
         times = times+1

      END DO RK

!
! End of Runge-Kutta

! Computes the benchmark

      IF (bench.gt.0) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu1)
         CALL GTStop(ihomp1)
         CALL GTStop(ihwtm1)
         inquire( file='benchmark.txt', exist=bbenchexist )
         IF (myrank.eq.0) THEN
            OPEN(1,file='benchmark.txt',position='append')
#if defined(DEF_GHOST_CUDA_)
            IF ( .NOT. bbenchexist ) THEN
               WRITE(1,*) &
	       '# nx ny nz nsteps nprocs nth nstrm TCPU TOMP TWTIME TFFT TTRA TCOM TMEM TASS TTOT'
            ENDIF
            WRITE(1,*) nx,ny,nz,(step-ini+1),nprocs,nth, &
                       nstreams                        , &
                       GTGetTime(ihcpu1)/(step-ini+1)  , &
                       GTGetTime(ihomp1)/(step-ini+1)  , &
                       GTGetTime(ihwtm1)/(step-ini+1)  , &
                       ffttime/(step-ini+1), tratime/(step-ini+1), &
                       comtime/(step-ini+1), memtime/(step-ini+1), &
                       asstime/(step-ini+1), tottime/(step-ini+1)
            WRITE(*,*) 'wtime=', GTGetTime(ihwtm1)/(step-ini+1),   &
	               ' fft=', ffttime/(step-ini+1),    &
		       ' transp=',tratime/(step-ini+1),  &
		       ' comm=',comtime/(step-ini+1),    &
		       ' mem=', memtime/(step-ini+1),    &
		       ' ttot=',tottime/(step-ini+1)
#else
            IF ( .NOT. bbenchexist ) THEN
               WRITE(1,*) &
	       '# nx ny nz nsteps nprocs nth TCPU TOMP TWTIME TFFT TTRA TCOM TTOT'
            ENDIF
            WRITE(1,*) nx,ny,nz,(step-ini+1),nprocs,nth, &
                       GTGetTime(ihcpu1)/(step-ini+1),   &
                       GTGetTime(ihomp1)/(step-ini+1),   &
                       GTGetTime(ihwtm1)/(step-ini+1),   &
                       ffttime/(step-ini+1), tratime/(step-ini+1), &
                       comtime/(step-ini+1), tottime/(step-ini+1)
            WRITE(*,*) 'wtime=', GTGetTime(ihwtm1)/(step-ini+1),   &
	               ' fft=', ffttime/(step-ini+1),    &
		       ' transp=',tratime/(step-ini+1),  &
		       ' comm=',comtime/(step-ini+1),    &
		       ' mem=',0.0, ' ttot=',tottime/(step-ini+1)
#endif
            IF (bench.eq.2) THEN
               WRITE(1,*) 'FFTW: Create_plan = ',      &
                       GTGetTime(ihcpu2)/(step-ini+1), &
                       GTGetTime(ihomp2)/(step-ini+1), &
                       GTGetTime(ihwtm2)/(step-ini+1)
            ENDIF
            CLOSE(1)
         ENDIF
      ENDIF

!
! End of MAIN3D

      CALL GTFree(ihcpu1)
      CALL GTFree(ihomp1)
      CALL GTFree(ihwtm1)
      CALL GTFree(ihcpu2)
      CALL GTFree(ihomp2)
      CALL GTFree(ihwtm2)

      CALL MPI_FINALIZE(ierr)
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)

      CONTAINS

!=================================================================
! SUBROUTINE: hdcheck_ndg
!
! Consistency check for the conservation of energy,
! helicity, and null divergence of the velocity field.
! Adapted from hdcheck in pseudospec3D_hd.f90 to take
! GStateComp field arrays and a filename label, for use
! with ensemble members and difference fields.
!
! Output files contain:
! '<label>balance.txt':    time, <v^2>, <omega^2>, injection rate
! '<label>helicity.txt':   time, kinetic helicity
! '<label>divergence.txt': time, <(div.v)^2>   [if chk=1]
!
! Parameters:
!   vel  : velocity field (3 components)
!   frc  : forcing field  (3 components)
!   t    : number of time steps made
!   dt   : time step
!   hel  : =0 skips helicity; =1 computes kinetic helicity
!   chk  : =0 skips divergence check; =1 performs it
!   label: prefix string for output filenames
!=================================================================
      SUBROUTINE hdcheck_ndg(vel, frc, t, dt, hel, chk, label)

          USE fprecision
          USE commtypes
          USE grid
          USE mpivars
          USE gstate_mod
          USE pseudospec_fluid
!$        USE threads

          IMPLICIT NONE

          TYPE(GStateComp), INTENT(IN) :: vel(:), frc(:)
          REAL(KIND=GP),    INTENT(IN) :: dt
          INTEGER,          INTENT(IN) :: t, hel, chk
          CHARACTER(LEN=*), INTENT(IN) :: label

          COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2,c3
          DOUBLE PRECISION :: eng,ens,pot,khe
          DOUBLE PRECISION :: div,tmp
          REAL(KIND=GP)    :: tmq
          INTEGER          :: i,j,k

          div = 0.0D0
          tmp = 0.0D0
          tmq = 1.0_GP/ &
                (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2

!
! Computes the mean square value of the divergence
!
          IF (chk.eq.1) THEN

          CALL derivk3(vel(1)%ccomp,c1,1)
          CALL derivk3(vel(2)%ccomp,c2,2)
          CALL derivk3(vel(3)%ccomp,c3,3)
          IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:tmp)
             DO j = 1,ny
                DO k = 1,nz
                   tmp = tmp+abs(c1(k,j,1)+c2(k,j,1)+c3(k,j,1))**2*tmq
                END DO
             END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:tmp)
             DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:tmp)
                DO j = 1,ny
                   DO k = 1,nz
                      tmp = tmp+2*abs(c1(k,j,i)+c2(k,j,i)+c3(k,j,i))**2*tmq
                   END DO
                END DO
             END DO
          ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:tmp)
             DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:tmp)
                DO j = 1,ny
                   DO k = 1,nz
                      tmp = tmp+2*abs(c1(k,j,i)+c2(k,j,i)+c3(k,j,i))**2*tmq
                   END DO
                END DO
             END DO
          ENDIF
          CALL MPI_REDUCE(tmp,div,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                          MPI_COMM_WORLD,ierr)

          ENDIF

!
! Computes mean energy, enstrophy, and kinetic helicity
!
          CALL energy(vel(1)%ccomp,vel(2)%ccomp,vel(3)%ccomp,eng,1)
          CALL energy(vel(1)%ccomp,vel(2)%ccomp,vel(3)%ccomp,ens,0)
          IF (hel.eq.1) THEN
             CALL helicity(vel(1)%ccomp,vel(2)%ccomp,vel(3)%ccomp,khe)
          ENDIF

!
! Computes the energy injection rate
!
          CALL cross(vel(1)%ccomp,vel(2)%ccomp,vel(3)%ccomp, &
                     frc(1)%ccomp,frc(2)%ccomp,frc(3)%ccomp,pot,1)

!
! Writes results to labelled output files
!
          IF (myrank.eq.0) THEN
             OPEN(1,file='balance'//TRIM(label)//'.txt',position='append')
             WRITE(1,10) (t-1)*dt,eng,ens,pot
   10        FORMAT( E13.6,E26.18,E26.18,E26.18 )
             CLOSE(1)
             IF (hel.eq.1) THEN
                OPEN(1,file='helicity'//TRIM(label)//'.txt',position='append')
                WRITE(1,FMT='(E13.6,E26.18)') (t-1)*dt,khe
                CLOSE(1)
             ENDIF
             IF (chk.eq.1) THEN
                OPEN(1,file='divergence'//TRIM(label)//'.txt',position='append')
                WRITE(1,FMT='(E13.6,E26.18)') (t-1)*dt,div
                CLOSE(1)
             ENDIF
          ENDIF

      END SUBROUTINE hdcheck_ndg

!=================================================================
! SUBROUTINE: hdspectrum_ndg
!
! Computes and writes the kinetic energy spectrum for a field
! given as a GStateComp array, appending a label to the filename.
! Output file: 'kspectrum.<ext><label>.txt'
!
! Parameters:
!   vel  : velocity field (3 components)
!   label: suffix appended to the spectrum index in the filename
!=================================================================
      SUBROUTINE hdspectrum_ndg(vel, label)

          USE fprecision
          USE mpivars
          USE gstate_mod
          USE filefmt
          USE pseudospec_fluid

          IMPLICIT NONE

          TYPE(GStateComp), INTENT(IN) :: vel(:)
          CHARACTER(LEN=*), INTENT(IN) :: label

          CALL spectrum(vel(1)%ccomp, vel(2)%ccomp, vel(3)%ccomp, &
                        TRIM(ext)//TRIM(label), 1, 1)

      END SUBROUTINE hdspectrum_ndg

!=================================================================
! SUBROUTINE: replace_scales
!
! Replaces the spectral modes of 'dst' with those of 'src' for
! all wavenumber shells k <= kndg. Modes at shells k > kndg in
! 'dst' are left unchanged. Both states must have the same number
! of components and be defined on the same distributed grid.
!
! Parameters:
!   src  : source field state (reference simulation)
!   dst  : destination field state (nudged simulation member)
!   kndg : cutoff wavenumber shell (integer)
!=================================================================
      SUBROUTINE replace_scales(src, dst, kcut)

          USE fprecision
          USE mpivars
          USE grid
          USE kes
          USE boxsize
          USE gstate_mod

          TYPE(GStateComp), INTENT(IN)    :: src(:)
          TYPE(GStateComp), INTENT(INOUT) :: dst(:)
          INTEGER,          INTENT(IN)    :: kcut

          INTEGER :: ic, i, j, k

          DO ic = 1, SIZE(src)
              DO i = ista, iend
                  DO j = 1, ny
                      DO k = 1, nz
                          IF (int(sqrt(kk2(k,j,i))/Dkk+0.5_GP) .le. kcut) THEN
                              dst(ic)%ccomp(k,j,i) = src(ic)%ccomp(k,j,i)
                          ENDIF
                      END DO
                  END DO
              END DO
          END DO

      END SUBROUTINE replace_scales

      END PROGRAM MAIN
