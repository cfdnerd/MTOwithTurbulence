# Gray Zone Issue – Root Cause Investigation

## Problem

The design stays at intermediate densities (gray) instead of converging to crisp 0/1. Tuning filterR/filterR2 can trigger or reduce this, but the underlying cause is deeper.

## Root Cause 1: Implementation Flaw (Primary)

**The cascaded filter uses two sequential *linear* PDE filters.** Each linear Helmholtz filter is a diffusion (smoothing) operation. Applying two in sequence produces *more* total diffusion than a single filter.

```
  x → [PDE 1: blur] → xp1 → [PDE 2: blur] → xp → [Heaviside] → xh
```

- **Two linear blurs** ⇒ combined effect = stronger smoothing ⇒ **more gray**, not less.
- Gray-suppressing cascaded filters in the literature (Lazarov & Wang, Hägg & Wadbro) use a **nonlinear** second stage (e.g. fW-mean with p>1, morphological operators, or projection) that *sharpens* or *projects*, not another linear blur.
- Our implementation: both stages are linear PDEs. The second stage adds diffusion, so gray increases. This is a **formulation bug**: the cascaded filter as implemented is not designed for gray suppression.

**Conclusion:** The legacy single-stage filter (one PDE + Heaviside) produces **less gray** because it applies only one blur before projection.

## Root Cause 2: Parameter Effects (Secondary)

Filter coefficients: `b = 1/(filterR·len/3.464)²`, `bFilter2 = 1/(filterR2·len/3.464)²`

- Larger radius → smaller `b` → more diffusion → more gray
- Increasing filterR or filterR2 (both linear stages) worsens gray
- Heaviside only sharpens the transition; it cannot remove gray if xp already has a wide band of intermediate values

## Fix Options

1. **Legacy filter**: `useCascadedFilter no` — single-stage PDE + Heaviside, avoids double-blur.
2. **Cascaded with nonlinear stage**: `useCascadedFilter yes` + `filterStage2Nonlinear yes` — fW-mean power stage suppresses gray (recommended when using cascaded).

## Nonlinear Second Stage (fW-mean)

With `filterStage2Nonlinear yes`, stage 2 applies xp = (W·xp1^p)^{1/p} instead of a linear PDE. For p>1 this sharpens (reduces gray). Parameters:
- **filterPMean** (default 4): p in the power mean; higher p = stronger gray suppression (4–8 typical).
- **filterR2**: PDE kernel radius for the fW-mean stage (neighborhood size).

## Parameter Guidance

- **filterR2** keep moderate (1.2–1.6) for stage-2 neighborhood
- **filterR** keep moderate (≤ 2.5); large values increase stage-1 gray
- **filterPMean** 4–8 for gray suppression without excessive sharpening
- **movlim** 0.3 helps reduce MMA oscillations
