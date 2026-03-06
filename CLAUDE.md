# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository extends the [GHOST](https://github.com/pmininni/GHOST) pseudospectral turbulence code (located at `~/repos/GHOST`) with a **spectral nudging / scale synchronization** ensemble approach. The idea is to run a reference simulation alongside one or more "nudged" simulations that are periodically overwritten at large scales (wavenumbers ≤ `kndg`) to match the reference, while small scales evolve freely.

## Building

This repo replaces `~/repos/GHOST/src/main.fpp`. To build, copy the files into the GHOST source tree and use GHOST's CMake build system:

```bash
cp main.fpp ~/repos/GHOST/src/
cd ~/repos/GHOST
cmake -B build -DSOLVER=HD -DPRECISION=DOUBLE -DFFTP=fftp-3
cmake --build build
```

Key CMake options (set in `UserConfig.cmake` or passed via `-D`):
- `SOLVER`: e.g. `HD`, `MHD`, `BOUSS` — controls which PDE solver is compiled
- `PRECISION`: `SINGLE` or `DOUBLE`
- `FFTP`: `fftp-3` (FFTW3) or `fftp-mkl`
- `ARBSIZE`: enable non-cubic box support
- `P_HYBRID`: enable MPI+OpenMP hybrid parallelization

Library paths can be passed as: `cmake .. -DCMAKE_PREFIX_PATH=/opt/openmpi -DFFTW3_ROOT=/opt/fftw`

## Running

```bash
mpirun -np <N> ~/repos/GHOST/bin/GHOST
```

The code reads all parameters from `parameter.inp` at runtime (Fortran namelists). See `~/repos/GHOST/src/README` and `~/repos/GHOST/src/examples/` for parameter file structure and examples. Binary output directories are set via the `&status` namelist.

## Architecture

### GHOST background (relevant modules)

All modules live in `~/repos/GHOST/src/pseudo/pseudospec3D_mod.f90`:
- `MODULE kes`: spectral arrays `kk2(nz,ny,ista:iend)` (squared wavenumber magnitudes), `kx/ky/kz`
- `MODULE boxsize`: physical domain parameters `Lx,Ly,Lz,Dkk` (Fourier shell width)
- `MODULE grid`: `nx,ny,nz`
- `MODULE fft`: FFT plan management

MPI decomposition variables (`ista`,`iend`,`ksta`,`kend`,`nprocs`,`myrank`) come from `mpivars`. Precision kind `GP` comes from `fprecision`.

Field state type `GStateComp` (from `~/repos/GHOST/src/memmgt/gstate_mod.f90`) wraps a single complex spectral component array: `ccomp(nz, ny, ista:iend)`. A full field is `TYPE(GStateComp), ALLOCATABLE :: field(:)` — one element per PDE component (e.g. 3 for velocity).

Factory subroutines (`init_pdes_from_file`, `init_ic_from_file`, `init_forcing_from_file`, `build_stepper_from_file`) construct polymorphic objects from `parameter.inp`.

### This repository

**`main.fpp`** — Fork of `~/repos/GHOST/src/main.fpp`. Additions:
- `TYPE :: GNudgedSim`: encapsulates all state for one nudged ensemble member (`field`, `field_nxt`, `force`, `pde`, `stepper`, `iclist`, `forcemethod`, `num_components`). Shares the reference simulation's `workspace`.
- `TYPE(GNudgedSim), ALLOCATABLE :: ensemble(:)`: array of nudged members, currently allocated with 1 element.
- Time loop advances the reference simulation via `stepper%step`, then for each ensemble member calls `replace_scales` followed by `ensemble(ii)%pde%timestep`.

Additional subroutines for this extension live in a `CONTAINS` block at the end of `main.fpp`:
- `replace_scales(src, dst, kndg)`: copies spectral modes at shells ≤ `kndg` from `src` into `dst`. Shell index computed as `int(sqrt(kk2(k,j,i))/Dkk + 0.5)`, consistent with GHOST's spectrum routines.

### Known issues (not yet fixed)
- `GNudgedSim%field` is declared as `TYPE(GState)` but should be `TYPE(GStateComp)` to match `GState_alloc` and `replace_scales`.
- `ensemble(ii)%field_nxt` is computed by `pde%timestep` but never assigned back — result is discarded each step.
- `ensemble(ii)%force` is allocated but `force` (reference) is passed to `pde%timestep` instead.
- `init_allstates` is never called for ensemble members.
