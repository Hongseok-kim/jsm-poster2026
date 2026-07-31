# Clone-Based Weighting for Future-Dependent Event Definitions

Poster for **JSM 2026** (Boston, MA) — Hongseok Kim, CSL Behring.

This repository hosts the poster and its supplementary materials. If you
scanned the QR code on the poster, you landed on
[`SUPPLEMENTARY.md`](SUPPLEMENTARY.md) — start there.

## What this is about

Standard time-to-event analysis assumes that, at any fixed time *t*, one can
determine whether the event has occurred by *t*. Some applied event
definitions break this assumption: deciding whether the event occurred at
*t_j* requires data from *later* time points. We call these
**future-dependent event definitions**.

A common instance is a **confirmed diagnosis defined as two consecutive
positive tests**. A subject's first positive at *t₁* is only ratified as an
event after a second positive at *t₂* — and if no second test arrives
before censoring, the first positive is *ambiguous*: it may or may not be a
real event. Standard Kaplan–Meier handling of these ambiguous cases is
biased.

This poster proposes a **clone-and-weight correction** for the
Kaplan–Meier estimator that resolves the ambiguity using a per-individual
confirmation probability *pᵢ*, and shows in simulation that the corrected
estimator tracks the true survival curve across realistic variation in
testing schedule and false-positive rate.

## The poster

- [`poster.pdf`](poster.pdf) — compiled poster (48 × 36 in, 3-column landscape)

## Supplementary materials

- [**`SUPPLEMENTARY.md`**](SUPPLEMENTARY.md) — **start here**: full reference list, detailed simulation settings and scenarios, the complete simulation code, and a discussion of the residual false-positive bias seen in Figure 2 of the poster
- [`docs/true_p_i_derivation.md`](docs/true_p_i_derivation.md) — closed-form derivation of *pᵢ* under exponential onset, sensitivity = 1, false-positive rate *ε*

## Code

- [`scripts/sim_core.R`](scripts/sim_core.R) — all simulation functions and parameters (sourced by the driver script below)
- [`scripts/Simulation.R`](scripts/Simulation.R) — main figures: KM curves, bias-vs-ε, bias-vs-gap

## Reproducing the figures

```r
# from project root
setwd("/path/to/jsm-poster2026")
source("scripts/Simulation.R")        # writes figures/processed/sim_*.pdf
```

Requires R with `survival`, `data.table`, and base graphics.

## Contact

Hongseok Kim — `kim.hongseok24 [at] gmail.com`
