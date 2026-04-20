!=================================================================
! MODULE: synchro_mod
!
! Module for the spectral nudging / scale synchronization
! ensemble extension. Defines the GNudgedSim type and its
! diagnostic methods (global, spectra) as well as helper
! routines (hdcheck_ndg, spectrum_ndg, replace_scales).
!=================================================================
      MODULE synchro_mod

      USE fprecision
      USE gstate_mod
      USE equationbase_mod
      USE gstepperbase_mod
      USE ic_factory
      IMPLICIT NONE

! Type encapsulating all state for a single nudged simulation member
      TYPE :: GNudgedSim
         TYPE(GStateComp),    ALLOCATABLE :: field(:),field_nxt(:)
         CLASS(EquationBase), ALLOCATABLE :: pde
         CLASS(GStepperBase), ALLOCATABLE :: stepper
         CLASS(icChain),      ALLOCATABLE :: iclist(:)
         INTEGER                          :: num_components
      CONTAINS
         PROCEDURE :: global  => gnudgedsim_global
         PROCEDURE :: spectra => gnudgedsim_spectra
      END TYPE GNudgedSim

      CONTAINS

!=================================================================
! SUBROUTINE: gnudgedsim_global
!
! Type-bound global diagnostic for a nudged ensemble member.
! Delegates to pde%global for the ensemble field (which writes
! to the member's own todir_), then computes the difference
! field (ensemble - reference) using workspace temporaries and
! writes its global quantities to 'difference.txt' in the same
! directory via hdcheck_ndg.
!=================================================================
      SUBROUTINE gnudgedsim_global(this, ref_field, force, t)

          USE status

          IMPLICIT NONE

          CLASS(GNudgedSim), INTENT(IN) :: this
          TYPE (GStateComp), INTENT(IN) :: ref_field(:), force(:)
          INTEGER,           INTENT(IN) :: t

          COMPLEX(KIND=GP), POINTER :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
          LOGICAL :: bret

          CALL this%pde%global(this%field_nxt, force, t)

          CALL this%pde%workspace_%get_complex_tmp(vx, bret)
          CALL this%pde%workspace_%get_complex_tmp(vy, bret)
          CALL this%pde%workspace_%get_complex_tmp(vz, bret)

          vx = this%field_nxt(1)%ccomp - ref_field(1)%ccomp
          vy = this%field_nxt(2)%ccomp - ref_field(2)%ccomp
          vz = this%field_nxt(3)%ccomp - ref_field(3)%ccomp

          CALL hdcheck_ndg(vx, vy, vz, &
                           force(1)%ccomp, force(2)%ccomp, force(3)%ccomp, &
                           t, dt, 0, 0, this%pde%todir_)

          CALL this%pde%workspace_%free_complex_tmp(vz)
          CALL this%pde%workspace_%free_complex_tmp(vy)
          CALL this%pde%workspace_%free_complex_tmp(vx)

      END SUBROUTINE gnudgedsim_global

!=================================================================
! SUBROUTINE: gnudgedsim_spectra
!
! Type-bound spectral diagnostic for a nudged ensemble member.
! Delegates to pde%spectra for the ensemble field (full output
! including energy transfer), then computes the difference
! spectrum (ensemble - reference) and writes it to
! 'dspectrum.XXX.txt' in the member's todir_.
!=================================================================
      SUBROUTINE gnudgedsim_spectra(this, ref_field)

          USE filefmt

          IMPLICIT NONE

          CLASS(GNudgedSim), INTENT(IN) :: this
          TYPE (GStateComp), INTENT(IN) :: ref_field(:)

          COMPLEX(KIND=GP), POINTER :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
          LOGICAL :: bret

          CALL this%pde%spectra(this%field_nxt)

          CALL this%pde%workspace_%get_complex_tmp(vx, bret)
          CALL this%pde%workspace_%get_complex_tmp(vy, bret)
          CALL this%pde%workspace_%get_complex_tmp(vz, bret)

          vx = this%field_nxt(1)%ccomp - ref_field(1)%ccomp
          vy = this%field_nxt(2)%ccomp - ref_field(2)%ccomp
          vz = this%field_nxt(3)%ccomp - ref_field(3)%ccomp

          CALL spectrum_ndg(vx, vy, vz, this%pde%todir_, ext)

          CALL this%pde%workspace_%free_complex_tmp(vz)
          CALL this%pde%workspace_%free_complex_tmp(vy)
          CALL this%pde%workspace_%free_complex_tmp(vx)

      END SUBROUTINE gnudgedsim_spectra

!=================================================================
! SUBROUTINE: spectrum_ndg
!
! Computes the kinetic energy spectrum and writes it to
! 'dspectrum.XXX.txt'. Same normalization as spectrum():
! E = sum[E(k).Dkk].
!=================================================================
      SUBROUTINE spectrum_ndg(a, b, c, path, nmb)

          USE kes
          USE grid
          USE mpivars
          USE filefmt
          USE boxsize
          USE pseudospec_hd

          IMPLICIT NONE

          COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
          CHARACTER(LEN=*), INTENT(IN) :: path, nmb

          DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Ek, Hk
          INTEGER :: i

          CALL spectrumc(a, b, c, 1, 0, Ek, Hk)

          IF (myrank.eq.0) THEN
             OPEN(1,file=TRIM(path)//'/dspectrum.'//nmb//'.txt')
             DO i = 1, nmax/2+1
                WRITE(1,FMT='(E13.6,E23.15)') Dkk*i, 0.5_GP*Ek(i)/Dkk
             END DO
             CLOSE(1)
          ENDIF

      END SUBROUTINE spectrum_ndg

!=================================================================
! SUBROUTINE: hdcheck_ndg
!
! Computes and writes global quantities (energy, enstrophy,
! injection rate) for a velocity field given as raw spectral
! arrays, writing to 'difference.txt' in the specified directory.
!
! Output files (written to todir):
!   'difference.txt'  : time, <v^2>, <omega^2>, injection
!   'helicity.txt'    : time, kinetic helicity     [if hel=1]
!   'divergence.txt'  : time, <(div.v)^2>          [if chk=1]
!=================================================================
      SUBROUTINE hdcheck_ndg(vx, vy, vz, fx, fy, fz, t, dt, hel, chk, todir)

          USE commtypes
          USE grid
          USE mpivars
          USE pseudospec_fluid
!$        USE threads

          IMPLICIT NONE

          COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend), INTENT(INOUT) :: vx,vy,vz
          COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend), INTENT(IN)    :: fx,fy,fz
          REAL(KIND=GP),    INTENT(IN) :: dt
          INTEGER,          INTENT(IN) :: t, hel, chk
          CHARACTER(LEN=*), INTENT(IN) :: todir

          COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend) :: c1,c2,c3
          DOUBLE PRECISION :: eng,ens,pot,khe
          DOUBLE PRECISION :: div,tmp
          REAL(KIND=GP)    :: tmq
          INTEGER          :: i,j,k

          div = 0.0D0
          tmp = 0.0D0
          tmq = 1.0_GP/ &
                (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2

          IF (chk.eq.1) THEN
             CALL derivk3(vx,c1,1)
             CALL derivk3(vy,c2,2)
             CALL derivk3(vz,c3,3)
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

          CALL energy(vx,vy,vz,eng,1)
          CALL energy(vx,vy,vz,ens,0)
          IF (hel.eq.1) THEN
             CALL helicity(vx,vy,vz,khe)
          ENDIF
          CALL cross(vx,vy,vz,fx,fy,fz,pot,1)

          IF (myrank.eq.0) THEN
             OPEN(1,file=TRIM(todir)//'/difference.txt',position='append')
             WRITE(1,10) (t-1)*dt,eng,ens,pot
   10        FORMAT( E13.6,E26.18,E26.18,E26.18 )
             CLOSE(1)
             IF (hel.eq.1) THEN
                OPEN(1,file=TRIM(todir)//'/helicity.txt',position='append')
                WRITE(1,FMT='(E13.6,E26.18)') (t-1)*dt,khe
                CLOSE(1)
             ENDIF
             IF (chk.eq.1) THEN
                OPEN(1,file=TRIM(todir)//'/divergence.txt',position='append')
                WRITE(1,FMT='(E13.6,E26.18)') (t-1)*dt,div
                CLOSE(1)
             ENDIF
          ENDIF

      END SUBROUTINE hdcheck_ndg

!=================================================================
! SUBROUTINE: replace_scales
!
! Replaces the spectral modes of 'dst' with those of 'src' for
! all wavenumber shells k <= kcut. Modes at shells k > kcut in
! 'dst' are left unchanged.
!=================================================================
      SUBROUTINE replace_scales(src, dst, kcut)

          USE mpivars
          USE grid
          USE kes
          USE boxsize

          IMPLICIT NONE

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

      END MODULE synchro_mod

!=================================================================
      PROGRAM MAIN
!=================================================================
! GHOST code: Geophysical High Order Suite for Turbulence
!
! Numerically integrates several fluid dynamics equations
! in 2 and 3 dimensions with periodic boundary conditions
! and external forcing. A pseudo-spectral method is used to
! compute spatial derivatives, while different time steppings
! can be used to evolve the system in the time domain.
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

! Modules
      USE commtypes
      USE filefmt
      USE iovar
      USE fft
      USE threads
      USE offloading
      USE boxsize
      USE status
      USE pstatus
      USE gtimer
      USE gbench
      USE ic_factory
      USE force_factory
      USE equation_factory
      USE particle_factory
      USE icp_factory
      USE stepper_factory
      USE synchro_mod
      IMPLICIT NONE

! Arrays for the field and particle states, workspace, I/O, and solver classes
      TYPE   (GStateComp), ALLOCATABLE :: field(:),field_nxt(:),force (:)
      TYPE  (GPStateComp), ALLOCATABLE :: part (:),part_nxt (:)
      TYPE   (GWorkspace)              :: workspace
      TYPE       (ioplan)              :: planio
      CLASS(EquationBase), ALLOCATABLE :: fluid
      CLASS(ParticleBase), ALLOCATABLE :: particle
      CLASS(GStepperBase), ALLOCATABLE :: stepper
      CLASS     (icChain), ALLOCATABLE :: iclist(:)
      CLASS  (forceChain), ALLOCATABLE :: forcemethod(:)
      CLASS    (icpChain), ALLOCATABLE :: icplist(:)

! Synchro ------------------------------------------
      TYPE(GNudgedSim), ALLOCATABLE :: ensemble(:)

! Auxiliary variables
      REAL(KIND=GP) :: time
      INTEGER       :: t, num_components

! Synchro ------------------------------------------
      INTEGER :: num_realizations, kndg
      INTEGER :: ir
      CHARACTER(LEN=8)   :: outlabel
      NAMELIST / ensembleparams / num_realizations, kndg
! Synchro ------------------------------------------

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

! Initialization of fluid and particles integration parameters
      CALL status_init ('parameter.inp')
      CALL pstatus_init('parameter.inp')

! Synchro ------------------------------------------
      num_realizations = 1
      kndg      = 1
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form='formatted')
         READ(1,NML=ensembleparams)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(num_realizations,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kndg     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! Synchro ------------------------------------------

! I/O initialization
      CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,(/nx,ny,nz/),ksta,kend,planio)

! Now we can initialize the PDE methods and read the particles status
      fluid        = init_pdes_from_file('parameter.inp')
      if (dopart) particle = init_particles_from_file('parameter.inp')
      CALL workspace%initialize_pool(NUMTMPREAL,NUMTMPCOMP,NUMTMPPART)
      CALL fluid%Solver_ctor('parameter.inp',workspace,planio)
      num_components = fluid%state_size()
      CALL GState_alloc(field    , num_components)
      CALL GState_alloc(field_nxt, num_components)
      CALL GState_alloc(force    , num_components)
      iclist       = init_ic_from_file(     'parameter.inp')
      forcemethod  = init_forcing_from_file('parameter.inp',workspace)

! Synchro ------------------------------------------
      ALLOCATE(ensemble(num_realizations))
      DO ir = 1, SIZE(ensemble)
         ! BLOCK and MOVE_ALLOC are needed becauce the base class
         ! does not initilize the pointers to null ("=> null()" in declaration)
         ! This is a workaround
         ! Maybe I can delete them now
         WRITE(outlabel,'(I3.3)') ir
         BLOCK
            CLASS(EquationBase), ALLOCATABLE :: tmp_pde
            tmp_pde = init_pdes_from_file('parameter'//TRIM(outlabel)//'.inp')
            CALL MOVE_ALLOC(tmp_pde, ensemble(ir)%pde)
         END BLOCK
         CALL ensemble(ir)%pde%Solver_ctor('parameter'//TRIM(outlabel)//'.inp', &
                                           workspace, planio)
         ensemble(ir)%num_components = ensemble(ir)%pde%state_size()
         CALL GState_alloc(ensemble(ir)%field    , ensemble(ir)%num_components)
         CALL GState_alloc(ensemble(ir)%field_nxt, ensemble(ir)%num_components)
         BLOCK
            CLASS(icChain),    ALLOCATABLE :: tmp_ic(:)
            tmp_ic  = init_ic_from_file('parameter'//TRIM(outlabel)//'.inp')
            CALL MOVE_ALLOC(tmp_ic,  ensemble(ir)%iclist)
         END BLOCK
      END DO
! Synchro ------------------------------------------

! Initialization of the numerical domain
      CALL box_init('parameter.inp')

! Initializes the FFT library. This must be done at this
! stage as it requires status and benchmark initialization.
! Use FFTW_ESTIMATE or FFTW_MEASURE in short runs
! Use FFTW_PATIENT or FFTW_EXHAUSTIVE in long runs
      nth = 1
!$    nth = omp_get_max_threads()
#if !defined(DEF_GHOST_CUDA_)
!$    CALL fftp3d_init_threads(ierr)
#endif
      CALL GTBenchInit(bench,ihcpu1,ihomp1,ihwtm1,ihcpu2,ihomp2,ihwtm2)
      IF (bench.eq.2) THEN
         CALL GTStart(ihcpu2); CALL GTStart(ihomp2); CALL GTStart(ihwtm2)
      ENDIF
      CALL fftp3d_create_plan(planrc,(/nx,ny,nz/),FFTW_REAL_TO_COMPLEX, &
                             FFTW_ESTIMATE)
      CALL fftp3d_create_plan(plancr,(/nx,ny,nz/),FFTW_COMPLEX_TO_REAL, &
                             FFTW_ESTIMATE)
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu2);  CALL GTStop(ihomp2);  CALL GTStop(ihwtm2)
      ENDIF

! Initial states
 IC : IF (stat.eq.0) THEN                 ! If stat=0 we start a new run
        ini  = 1
        sind = 0                          ! index for the spectrum
        tind = 0                          ! index for the binaries
        pind = 0                          ! index for the particles
        timet = tstep
        timep = pstep
        timec = cstep
        times = sstep
        timef = 0
      ELSE                    ! If stat.ne.0 a previous run is continued
        ini = int((stat-1)*tstep) + 1
        tind = int(stat)
        sind = int(real(ini,kind=GP)/real(sstep,kind=GP)+1)
        pind = int((stat-1)*lgmult+1)
        timet = 0
        timep = 0
        timec = int(modulo(float(ini-1),float(cstep)))
        times = int(modulo(float(ini-1),float(sstep)))
        timef = int(modulo(float(ini-1),float(fstep)))
      ENDIF IC
      CALL init_allstates(iclist,fluid,field)
      field_nxt = field  ! We update nxt to work with I/O and all steppers
      CALL init_forcing(forcemethod,fluid,force)
      if (dopart) then
         CALL init_allpstates(icplist,fluid,field,particle,part)
	 if (size(part(1)%rcomp) .ne. size(part_nxt(1)%rcomp)) then
	   call GPState_resize(part_nxt,particle%partbuff_) ! We resize part_nxt
	 endif
	 part_nxt = part ! We also update part_nxt
     endif

! Synchro ------------------------------------------
! Initialize states and stepper
      DO ir = 1, SIZE(ensemble)
      WRITE(outlabel,'(I3.3)') ir
      CALL init_allstates(ensemble(ir)%iclist, ensemble(ir)%pde, ensemble(ir)%field)
      BLOCK
         CLASS(GStepperBase), ALLOCATABLE :: tmp_stp
         tmp_stp = build_stepper_from_file('parameter'//TRIM(outlabel)//'.inp', &
                                           workspace, &
                                           ensemble(ir)%pde)
         CALL MOVE_ALLOC(tmp_stp, ensemble(ir)%stepper)
      END BLOCK
      END DO

! Initial replace scales
      DO ir = 1, SIZE(ensemble)
         CALL replace_scales(field, ensemble(ir)%field, kndg)
         ensemble(ir)%field_nxt = ensemble(ir)%field
      END DO
! Synchro ------------------------------------------

! Sets up the time stepper
      stepper = build_stepper_from_file('parameter.inp',workspace,fluid)

! Time integration scheme starts here.
! If we are doing a benchmark, we measure cputime before
! starting. We also re-inititialize the fftp timers.
      IF (bench.eq.1) THEN
          ffttime  = 0.D00; tratime  = 0.0D0; comtime  = 0.D00; tottime  = 0.0D0
          CALL GTStart(ihcpu1); CALL GTStart(ihomp1); CALL GTStart(ihwtm1)
      ENDIF

RK :  DO t = ini,step
         time = (t-1)*dt
! Every 'tstep' steps, stores the fields in binary files
         IF ((timet.eq.tstep).and.(bench.eq.0)) THEN
            timet = 0
            tind = tind+1
            CALL fluid%write_states(field_nxt, planio)
            DO ir = 1, SIZE(ensemble)
               CALL ensemble(ir)%pde%write_states(ensemble(ir)%field_nxt, planio)
            END DO
         ENDIF

! Every 'cstep' steps writes global quantities
         IF ((timec.eq.cstep).and.(bench.eq.0)) THEN
            timec = 0
     	      CALL fluid%global(field_nxt, force, t)
            DO ir = 1, SIZE(ensemble)
               CALL ensemble(ir)%global(field_nxt, force, t)
            END DO
         ENDIF

! Every 'sstep' steps writes spectra
         IF ((times.eq.sstep).and.(bench.eq.0)) THEN
            times = 0
            sind = sind+1
            CALL fluid%spectra(field_nxt)
            DO ir = 1, SIZE(ensemble)
               CALL ensemble(ir)%spectra(field_nxt)
            END DO
         ENDIF

! Time evolution
         CALL update_forcing(forcemethod,fluid,force)
         field = field_nxt
         CALL stepper%gstep(time, field, force, dt, field_nxt)
! Synchro ----------------------------------------------------------
! step ensemble members and replace large scales
         DO ir = 1, SIZE(ensemble)
            ensemble(ir)%field = ensemble(ir)%field_nxt
            CALL ensemble(ir)%stepper%gstep(time, ensemble(ir)%field, &
                                            force, dt, ensemble(ir)%field_nxt)
            CALL replace_scales(field_nxt, ensemble(ir)%field_nxt, kndg)
         END DO
! Synchro ----------------------------------------------------------
         timet = timet+1; timep = timep+1; timec = timec+1; times = times+1
      END DO RK

! Finishes and writes the benchmark results
      IF (bench.gt.0) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu1); CALL GTStop(ihomp1); CALL GTStop(ihwtm1)
      ENDIF
      CALL GTBenchReport(bench,myrank,ini,step,ihcpu1,ihomp1,ihwtm1,             &
               ihcpu2,ihomp2,ihwtm2)
      IF (dopart) THEN
        rbal = rbal + GetLoadBal(particle) ! Get load balancing
        CALL GTBenchReportParticles(bench, myrank, ini, step, maxparts, rbal,    &
             [ GetTime(particle,GPTIME_STEP),   GetTime(particle,GPTIME_COMM),   &
	       GetTime(particle,GPTIME_SPLINE), GetTime(particle,GPTIME_TRANSP), &
	       GetTime(particle,GPTIME_DATAEX), GetTime(particle,GPTIME_INTERP), &
	       GetTime(particle,GPTIME_PUPDATE) ], 'gpbenchmark.txt')
      ENDIF

! End of main
      CALL GTFree(ihcpu1)
      CALL GTFree(ihomp1)
      CALL GTFree(ihwtm1)
      CALL GTFree(ihcpu2)
      CALL GTFree(ihomp2)
      CALL GTFree(ihwtm2)
      CALL MPI_FINALIZE(ierr)
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)

      END PROGRAM MAIN
