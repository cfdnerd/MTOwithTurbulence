# MTO_ThermalFluid Case Setup

## Symmetric boundary condition

This case is configured for a **symmetric half-domain** geometry:

- **Right boundary** (`right` patch): symmetry plane — the mesh represents one half of the full domain; the other half is the mirror image.
- All fields (U, p, T, Ua, Ub, pa, pb, Tb, k, omega, ka, omegaa) use `type symmetry` on the `right` patch.

### Implications

1. **Power dissipation**: The cost function integrates over inlet and outlet. With symmetry, the computed power dissipation corresponds to the half-domain; PowerDiss0 and constraints should be set accordingly.
2. **Mesh**: Run `blockMesh` to generate the half-domain mesh. The `right` patch must remain a symmetry plane (no flow normal to it, zero tangential stress).

### Inlet / outlet

- Inlet velocity: 1 m/s
- Inlet port diameter: 2 mm (slot width in 2D)
