# Cascaded PDE Filter for Turbulent Channel Flow Topology Optimization

## Overview

The cascaded filter supports two modes for stage 2:
- **Linear PDE** (filterStage2Nonlinear no): Second Helmholtz blur — can increase gray.
- **Nonlinear fW-mean** (filterStage2Nonlinear yes): Power-mean stage — suppresses gray (recommended).

## Implementation

### Filter Pipeline

```
  x  →  Stage 1 (PDE)  →  xp1  →  Stage 2 (linear or fW-mean)  →  xp  →  Heaviside  →  xh
```

1. **Stage 1**: Helmholtz PDE with radius `filterR` — `laplacian(xp1) - b*xp1 + b*x = 0`
2. **Stage 2** (optional nonlinear):
   - Linear: `laplacian(xp) - bFilter2*xp + bFilter2*xp1 = 0`
   - Nonlinear (fW-mean): `xp = (W·xp1^p)^{1/p}` using PDE kernel W (radius filterR2)
3. **Heaviside**: Volume-preserving projection with `del = min(0.4*opt, 200)`

### Parameters (optProperties)

| Parameter             | Default | Description                                      |
|-----------------------|---------|--------------------------------------------------|
| filterR               | 2.2     | Stage 1 filter radius                            |
| filterR2              | 1.6     | Stage 2 PDE kernel radius                        |
| filterStage2Nonlinear | yes     | Use fW-mean (gray suppression) vs linear PDE     |
| filterPMean           | 4       | fW-mean exponent p (p>1; 4–8 typical for gray suppression) |

### Sensitivity Chain Rule

The adjoint sensitivities pass through:

1. Heaviside derivative: `drho = d(xh)/d(xp)`
2. Stage 2 PDE adjoint: `dJ/dxp → dJ/dxp1`
3. Stage 1 PDE adjoint: `dJ/dxp1 → dJ/dx`

## Legacy Filter (Recommended for Gray Suppression)

The cascaded filter uses two linear PDE blurs, which increases gray. Use the legacy single-stage filter when gray is an issue:

**Set in optProperties:** `useCascadedFilter no;`  (default if omitted)

To switch back manually (instead of using the switch):
1. Replace `#include "filter_chainrule.H"` with `#include "filter_chainrule_legacy.H"` in `sensitivity.H`
2. Replace `#include "filter_x.H"` with `#include "filter_x_legacy.H"` in `sensitivity.H`
3. Remove or comment out `filterR2` and `bFilter2` usage if desired (legacy ignores them)

## References

- Lazarov, B.S., Wang, F. (2016). "Nonlinear filters in topology optimization: existence of solutions and efficient implementation for minimum compliance problems." *Struct Multidisc Optim* 54, 15–32.
- Guest et al. (2006). "Morphology-based black and white filters for topology optimization." *Struct Multidisc Optim* 33, 401–424.
