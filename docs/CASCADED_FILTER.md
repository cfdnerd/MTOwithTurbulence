# Cascaded PDE Filter for Turbulent Channel Flow Topology Optimization

## Overview

The cascaded density filter (Lazarov & Wang, 2016) has been implemented to improve gray suppression and 0/1 design convergence in turbulent internal channel flows. The previous single-stage PDE filter is preserved in `filter_x_legacy.H` and `filter_chainrule_legacy.H`.

## Implementation

### Filter Pipeline

```
  x (design)  →  Stage 1 PDE  →  xp1  →  Stage 2 PDE  →  xp  →  Heaviside  →  xh
```

1. **Stage 1**: Helmholtz PDE filter with radius `filterR`  
   - Solves: `laplacian(xp1) - b*xp1 + b*x = 0`
2. **Stage 2**: Second Helmholtz PDE filter with radius `filterR2`  
   - Solves: `laplacian(xp) - bFilter2*xp + bFilter2*xp1 = 0`
3. **Heaviside**: Volume-preserving projection with `del = min(0.4*opt, 200)`

### Parameters (optProperties)

| Parameter  | Default | Description                                        |
|------------|---------|----------------------------------------------------|
| filterR    | 2.2     | Stage 1 filter radius (design smoothing)           |
| filterR2   | 1.6     | Stage 2 filter radius (smaller = sharper features) |

Using `filterR2 < filterR` produces stronger gray suppression while maintaining length-scale control.

### Sensitivity Chain Rule

The adjoint sensitivities pass through:

1. Heaviside derivative: `drho = d(xh)/d(xp)`
2. Stage 2 PDE adjoint: `dJ/dxp → dJ/dxp1`
3. Stage 1 PDE adjoint: `dJ/dxp1 → dJ/dx`

## Legacy Filter

To revert to the original single-stage filter:

1. In `sensitivity.H`, replace `#include "filter_chainrule.H"` with `#include "filter_chainrule_legacy.H"`
2. Replace `#include "filter_x.H"` with `#include "filter_x_legacy.H"`
3. Remove or comment out the `filterR2` and `bFilter2` usage in `opt_initialization.H` and `createFields.H` (or keep them; the legacy filter ignores them)
4. Remove the `xp1` field if desired (legacy does not use it)

## References

- Lazarov, B.S., Wang, F. (2016). "Nonlinear filters in topology optimization: existence of solutions and efficient implementation for minimum compliance problems." *Struct Multidisc Optim* 54, 15–32.
- Guest et al. (2006). "Morphology-based black and white filters for topology optimization." *Struct Multidisc Optim* 33, 401–424.
