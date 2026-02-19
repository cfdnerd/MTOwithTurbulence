# MTO_ThermalFluid — Build for OpenFOAM 6

## Environment

- **OpenFOAM:** 6
- **MPI:** OpenMPI 1.10.7 (must match the MPI used to build OpenFOAM and PETSc)
- **PETSc:** Required; configure with OpenMPI 1.10.7

## Build Steps

1. Source OpenFOAM 6:
   ```bash
   source /opt/openfoam6/etc/bashrc
   # or wherever OpenFOAM 6 is installed
   ```

2. Ensure `WM_PROJECT_USER_DIR` is set (e.g. `$HOME/OpenFOAM/$USER-6`).

3. Build MTO_ThermalFluid:
   ```bash
   cd $MTO_SOURCE/MTO_ThermalFluid
   wmake
   ```

4. Executable will be placed in `$FOAM_USER_APPBIN/MTO_ThermalFluid`.

## Make/options Notes

- Uses `$(LIB_SRC)` for OpenFOAM library paths (set by wmake).
- PETSc requires `PETSC_DIR` and `PETSC_ARCH` set in the environment before building.

## MPI Compatibility

If PETSc is built with OpenMPI 1.10.7, build OpenFOAM 6 with the same OpenMPI to avoid ABI conflicts. Use `mpirun` from that OpenMPI for parallel runs.
