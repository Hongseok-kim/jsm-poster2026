# Supplementary Materials

**Clone-Based Weighting for Future-Dependent Event Definitions**
Poster — JSM 2026 (Boston, MA). Hongseok Kim, CSL Behring.
[LinkedIn](https://www.linkedin.com/in/hongseok-kim1/)

This document collects the material that does not fit on the poster: the full
reference list, a detailed description of the simulation settings and scenarios,
the complete simulation code, and a derivation of the proposed estimator from
the standard Kaplan–Meier formula. It closes with a discussion of where the
confirmation probability *pᵢ* comes from, and of the residual false-positive
bias visible in Figure 2 of the poster. Inference for the estimator is still
under development, as §4 notes.

**Contents**

1. [References](#1-references)
2. [Simulation settings and scenarios](#2-simulation-settings-and-scenarios)
3. [Simulation code](#3-simulation-code)
4. [Inference for the proposed estimator](#4-inference-for-the-proposed-estimator)
5. [Estimating *pᵢ*](#5-estimating-pᵢ)
6. [Residual false-positive bias and its correction](#6-residual-false-positive-bias-and-its-correction)

Related notes already in the repository:

- [`docs/true_p_i_derivation.md`](docs/true_p_i_derivation.md) — closed-form derivation of the oracle weight *pᵢ*.

---

## 1. References

The poster cites the following works, listed here in order of first appearance
on the poster. Numbers are for this document only; the poster itself uses
author–year citations.

**Foundations of time-to-event analysis**

1. Andersen PK, Borgan Ø, Gill RD, Keiding N. *Statistical Models Based on
   Counting Processes.* Springer Series in Statistics. New York: Springer; 1993.

**Confirmatory / repeated-testing event definitions in practice**

2. Saracino A, et al. Increased risk of virologic failure to the first
   antiretroviral regimen in HIV-infected migrants compared to natives: data
   from the ICONA cohort. *Clinical Microbiology and Infection.*
   2016;22(3):288.e1–288.e9.

3. Phillips PPJ, et al. Limited role of culture conversion for decision-making
   in individual patient care and for advancing novel regimens to confirmatory
   clinical trials. *BMC Medicine.* 2016;14(1):19.

4. Heerspink HJL, et al. Dapagliflozin in patients with chronic kidney disease.
   *New England Journal of Medicine.* 2020;383(15):1436–1446.

**Confirmation reduces false positives**

5. Selvin E, et al. Prognostic implications of single-sample confirmatory
   testing for undiagnosed diabetes: a prospective cohort study. *Annals of
   Internal Medicine.* 2018;169(3):156–164.

6. Moses MW, et al. Serial testing for latent tuberculosis using QuantiFERON-TB
   Gold In-Tube: a Markov model. *Scientific Reports.* 2016;6(1):30781.

---

## 2. Simulation settings and scenarios

The simulation compares four Kaplan–Meier estimators on synthetic cohorts of
repeated diagnostic testing, where the event of interest ("confirmed diagnosis")
is defined as **two consecutive positive tests**. The goal is to show that the
proposed clone-weighted estimator tracks the true onset survival curve across
realistic variation in the testing schedule and the false-positive rate.

### 2.1 Scope and assumptions

The whole construction is developed under the following deliberate scope choices
(not oversights).

- **Sensitivity = 1 (no false negatives) — the key assumption.** A truly
  positive person tests positive at every test from onset onward. This makes an
  observed negative *certainly* truly negative (`RE = 0 ⟹ R = 0`) and gives each
  true case an unbroken run of positives. Relaxing this to sensitivity < 1 is
  future work.
- **Test error is false positives only**, at rate
  ε = P(test + | truly negative) = 1 − specificity. Since specificity is
  routinely reported for diagnostic assays, ε is often known in practice.
- **False positives are independent across tests** within a person.

### 2.2 Data-generating model

For each individual:

- **True onset time** *T₀* ~ Exp(λ), so the true survival curve is
  *S(t) = e^(−λt)*. This is the target the estimators are judged against.
- **Test times** *O₁ < O₂ < … < O_K*. The first test is drawn as
  *O₁ ~ N(1, gap_sd)*. Each subsequent gap depends on the **observed** result of
  the previous test: shorter after a positive (closer monitoring), longer after
  a negative (routine interval). Scheduling follows what the clinician sees, so a
  false positive also triggers the shorter follow-up.
- **True status** at test *k*: *R_k = 1{T₀ ≤ O_k}*.
- **Observed result** at test *k*: equal to *R_k* when truly positive
  (sensitivity = 1); when truly negative, a false positive occurs independently
  with probability ε.

An event under the confirmed-diagnosis rule is the **first of two consecutive
observed positives**, dated at the first of the pair. A **censored
first-positive** is an individual with a first positive before the analysis
horizon whose confirming (next) test falls beyond the horizon — the ambiguous
case the proposed method resolves.

### 2.3 Fixed parameters

Set in [`scripts/sim_core.R`](scripts/sim_core.R) (`set.seed(2026)`):

| Parameter | Symbol | Value | Meaning |
|---|---|---|---|
| Sample size | *n* | 1000 | individuals per replicate |
| Onset rate | λ | 0.15 | *T₀ ~ Exp(λ)* |
| Tests per person | *K* | 5 | scheduled tests |
| Horizon / cutoff | — | 2 | administrative censoring & comparison time |
| Replicates | — | 500 | Monte-Carlo repetitions per setting |
| Gap SD | gap_sd | 0.1 | common SD of inter-test gaps |
| After-positive gap (baseline) | gap_pos | 0.40 | mean gap after an observed positive |
| After-negative gap (baseline) | gap_neg | 0.75 | mean gap after an observed negative |
| False-positive rate (baseline) | ε | 0.05 | 1 − specificity |

### 2.4 The four estimators

| Estimator | Event rule | Role |
|---|---|---|
| **True** | first positive on the *error-free* status *R* | gold standard (uses the unobservable truth) |
| **Single-positive** | first *observed* positive | naive; over-counts events (false positives) |
| **Double-positive** | first of two consecutive observed positives; censored first-positives left as non-events | naive; under-counts events (drops ambiguous cases) |
| **Clone-weighted** | double-positive, but each censored first-positive is split into an event clone (weight *pᵢ*) and a censored clone (weight 1 − *pᵢ*) | **proposed** |

The clone-weighted estimator uses the **oracle** weight *pᵢ* from the closed
form in [`docs/true_p_i_derivation.md`](docs/true_p_i_derivation.md), which
isolates the clone-weighting *mechanism* from the separate problem of
*estimating* *pᵢ* (see §5).

### 2.5 Scenarios

One parameter is varied at a time; the other is held at its baseline. Bias is
measured on the event rate at the cutoff, signed so that **positive = the
estimator overstates the event rate** and **negative = understates**:

> bias = (1 − S_est) − (1 − S_true) = S_true − S_est.

**Baseline KM curves** (`sim_main_km.pdf`) — all four curves at
ε = 0.05, gap_pos = 0.40. Illustrates the qualitative pattern: Single overstates,
Double understates, Clone tracks the truth.

**Scenario A — false-positive sweep** (`sim_bias_vs_eps.pdf`):
ε ∈ {0.02, 0.05, 0.10, 0.15}, gap_pos fixed at 0.40. Isolates the effect of
false-positive contamination.

**Scenario B — schedule sweep** (`sim_bias_vs_gap.pdf`):
after-positive gap ∈ {0.25, 0.40, 0.55, 0.75}, ε fixed at 0.05 (after-negative
gap held at 0.75). As the confirming test moves later, more genuine
first-positives become censored before they can be confirmed. At the top of the
grid the after-positive gap equals the after-negative gap, i.e. no closer
monitoring after a positive.

In Scenario A the clone-weighted estimator does not sit exactly at zero, and its
residual bias grows with ε. That residual is a property of the double-positive
event definition rather than of clone weighting, and it is discussed in §6.

### 2.6 Reproducing the figures

```r
# from the project root
setwd("/path/to/jsm-poster2026")
source("scripts/Simulation.R")        # writes figures/processed/sim_main_km.pdf,
                                       #        sim_bias_vs_eps.pdf, sim_bias_vs_gap.pdf
```

Requires R with the `survival` and `data.table` packages and base graphics.
Running from the project root ensures the figures land in the canonical
`figures/processed/` directory.

---

## 3. Simulation code

The code below is the canonical source in [`scripts/`](scripts/), reproduced here
so the supplement is self-contained. If it ever disagrees with the files in
`scripts/`, the files are authoritative.

### 3.1 `sim_core.R` — shared functions and parameters

```r
# ============================================================================
# Simulation -- shared core
#
# Functions, fixed parameters, and plotting helpers used by:
#   * Simulation.R       -- main figures (True / Single / Double / Clone)
#
# This file only DEFINES things and sets parameters; it runs no sweeps and
# writes no figures. Source it from a driver script:
#     source(if (file.exists("sim_core.R")) "sim_core.R" else "scripts/sim_core.R")
#
# ASSUMPTIONS (scope of this work):
#   * SENSITIVITY = 1 (no false negatives): a truly positive person tests
#     positive at every test from onset on. This is the key simplification --
#     it makes an observed negative certainly truly negative (RE=0 => R=0) and
#     gives each true case an unbroken run of positives. Extension to
#     sensitivity < 1 is future work.
#   * Test error is false positives only, at rate epsilon = 1 - specificity,
#     assumed KNOWN (commonly available as reported assay specificity).
#   * False positives are independent across tests within a person.
# ============================================================================

library(survival)
library(data.table)

set.seed(2026)

# -- Fixed parameters --------------------------------------------------------
n          <- 1000
lambda1    <- 0.15
n_tests    <- 5
cutoff     <- 2
n_reps     <- 500
gap_sd     <- 0.1

# Baselines (used when sweeping the other parameter)
gap_pos_baseline <- 0.4    # mean gap after an observed positive (closer monitoring)
gap_neg_baseline <- 0.75   # mean gap after an observed negative (routine interval)
eps_baseline     <- 0.05

# Sweep grids
eps_grid <- c(0.02, 0.05, 0.10, 0.15)
# gap sweep varies the after-positive gap up to gap_neg_baseline (= 0.75), where
# there is no longer any closer monitoring after a positive; includes the
# baseline 0.4 as a cross-check against the epsilon sweep's baseline point.
gap_grid <- c(0.25, 0.40, 0.55, 0.75)

# Output: save plots to figures/processed/ as PDFs?
SAVE_PLOTS <- TRUE
out_dir    <- file.path("figures", "processed")
if (SAVE_PLOTS) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -- Helpers -----------------------------------------------------------------

# ----------------------------------------------------------------------------
# gap_gen()
# Purpose: draw the time gap to the next test, which depends on the previous
#          result -- shorter after a positive (closer monitoring), longer
#          otherwise (routine monitoring). Defined before simulate_data(),
#          which calls it once per test from the second test onward.
# Args:
#   res_prev - integer 0/1 vector of the previous test's OBSERVED result, i.e.
#              the column RE[, k-1] passed in by simulate_data(), one entry per
#              individual. Scheduling follows what is actually observed, so a
#              false positive (RE = 1) also triggers the shorter mean_pos gap.
#              It is the branch driver: a 1 selects the shorter mean_pos gap, a
#              0 selects the longer mean_neg gap.
#   mean_pos - mean gap for individuals whose previous result was positive.
#   mean_neg - mean gap for everyone else.
#   sd       - common standard deviation of the gaps (default gap_sd).
# Returns: numeric vector of inter-test gaps, one per individual.
gap_gen <- function(res_prev, mean_pos, mean_neg, sd = gap_sd) {
  out       <- numeric(length(res_prev))
  pos       <- res_prev == 1L
  out[pos]  <- rnorm(sum(pos),  mean = mean_pos, sd = sd)
  out[!pos] <- rnorm(sum(!pos), mean = mean_neg, sd = sd)
  out
}

# ----------------------------------------------------------------------------
# simulate_data()
# Purpose: generate one simulated cohort of repeated diagnostic testing.
# Args:
#   n       - number of individuals.
#   lambda  - rate of the exponential true-onset time, T0 ~ Exp(lambda).
#   n_tests - number of scheduled tests per individual.
#   gap_pos - mean inter-test gap after a positive result (closer monitoring).
#   gap_neg - mean inter-test gap after a negative result (routine monitoring).
#   eps     - false-positive rate of a test when the truth is still negative.
# Returns: a list of three n x n_tests matrices -- O, R, RE (described below).
simulate_data <- function(n, lambda, n_tests, gap_pos, gap_neg, eps) {
  T0 <- rexp(n, rate = lambda)

  # Each matrix is n individuals (rows) x n_tests (columns):
  #   O  : test (observation) times -- when each test is taken. The first
  #        test ~ N(1, gap_sd); each later test is spaced from the previous one
  #        by gap_gen() keyed on the previous OBSERVED result RE (closer
  #        monitoring after an observed positive, false positives included).
  #   R  : true status (error-free). R[i, k] = 1 if the true onset T0
  #        has occurred by test time O[i, k], else 0. Sensitivity = 1 and
  #        there are no false positives, so R reflects the latent truth.
  #   RE : observed test results (with error). Equals R whenever truly
  #        positive; when truly negative, a false positive occurs with
  #        probability eps. RE is what an analyst actually sees with observed test.
  O  <- matrix(NA_real_,    nrow = n, ncol = n_tests)
  R  <- matrix(NA_integer_, nrow = n, ncol = n_tests)
  RE <- matrix(NA_integer_, nrow = n, ncol = n_tests)

  O[, 1]  <- rnorm(n, mean = 1, sd = gap_sd)
  R[, 1]  <- as.integer(T0 <= O[, 1])
  RE[, 1] <- ifelse(R[, 1] == 1L, 1L, rbinom(n, 1, eps))

  for (k in 2:n_tests) {
    O[, k]  <- O[, k - 1] + gap_gen(RE[, k - 1], gap_pos, gap_neg)
    R[, k]  <- as.integer(T0 <= O[, k])
    RE[, k] <- ifelse(R[, k] == 1L, 1L, rbinom(n, 1, eps))
  }

  list(O = O, R = R, RE = RE)
}

# ----------------------------------------------------------------------------
# compute_single()
# Purpose: turn a test-result matrix into a (time, event) survival outcome
#          under the SINGLE-positive rule -- the event happens at the first
#          positive test. Called twice: with R it yields the error-free "True"
#          outcome; with RE it yields the naive "Single-positive" outcome.
# Args:
#   O       - n x K matrix of test (observation) times.
#   results - n x K matrix of 0/1 test results (R for True, RE for Single).
#   cutoff  - administrative censoring time and analysis horizon.
# Returns: data.table with one row per individual --
#   T - observed time (first-positive time, or cutoff if censored).
#   E - event indicator (1 = event seen by cutoff, 0 = censored).
compute_single <- function(O, results, cutoff) {
  n_obs <- nrow(O); K <- ncol(O)
  T_out <- rep(NA_real_, n_obs); E_out <- rep(NA_integer_, n_obs)
  for (k in 1:K) {
    cens <- is.na(T_out) & (cutoff < O[, k])
    T_out[cens] <- cutoff;   E_out[cens] <- 0L
    evt  <- is.na(T_out) & (results[, k] == 1L)
    T_out[evt]  <- O[evt, k]; E_out[evt]  <- 1L
  }
  data.table(T = T_out, E = E_out)
}

# ----------------------------------------------------------------------------
# compute_double()
# Purpose: turn the observed results into a (time, event) outcome under the
#          TWO-consecutive-positives (confirmed) rule, with the event time set
#          to the FIRST of the two positives. Also records, for each individual,
#          the interval of a "censored first-positive" -- a first positive that
#          falls before the cutoff but whose confirming (next) test is beyond
#          the horizon, so it can never be confirmed in time. These are the
#          ambiguous cases that clone-weighting resolves.
# Args:
#   O      - n x K matrix of test (observation) times.
#   RE     - n x K matrix of observed (error-prone) 0/1 results.
#   cutoff - administrative censoring time and analysis horizon.
# Returns: a list with --
#   dt      - data.table(T, E) under the double-positive rule.
#   amb     - integer vector of length n; amb[j] is the interval (test index) of
#             individual j's censored first-positive, or NA if there is none.
#             Covers ALL intervals, including test 1 (previously only test 2 was
#             handled). Interval K is out of scope (no confirming test exists),
#             but does not arise for cutoff = 2 in this design.
#   evt_int - integer vector of length n; evt_int[j] is the interval (first
#             positive of the confirmed pair) at which individual j's event was
#             recorded, or NA if not a confirmed event.
compute_double <- function(O, RE, cutoff) {
  n_obs <- nrow(O); K <- ncol(O)
  T_out <- rep(NA_real_, n_obs); E_out <- rep(NA_integer_, n_obs)
  amb     <- rep(NA_integer_, n_obs)   # interval of a censored first-positive
  evt_int <- rep(NA_integer_, n_obs)   # interval of a confirmed event (else NA)

  # Interval of each individual's FIRST observed positive (0 if never positive).
  first_pos <- max.col(RE, ties.method = "first")
  first_pos[rowSums(RE) == 0L] <- 0L

  # Initial censoring: when the second test is beyond the horizon, a first
  # positive at test 1 can never be confirmed in time -> flag it (interval 1).
  cens <- cutoff < O[, 2]
  amb[cens & first_pos == 1L & O[, 1] < cutoff] <- 1L
  T_out[cens] <- cutoff; E_out[cens] <- 0L

  for (k in 2:K) {
    # Confirmed event: first two-in-a-row positives, dated at the first positive
    # (interval k-1).
    evt <- is.na(T_out) & RE[, k - 1] == 1L & RE[, k] == 1L
    T_out[evt] <- O[evt, k - 1]; E_out[evt] <- 1L
    evt_int[evt] <- k - 1L

    if (k < K) {
      # Censor those whose confirming test (k+1) is beyond the horizon; among
      # them, a first positive at test k (before the cutoff) is a censored
      # first-positive at interval k -> flag it for cloning.
      cens <- is.na(T_out) & (cutoff < O[, k + 1])
      amb[cens & first_pos == k & O[, k] < cutoff] <- k
      T_out[cens] <- cutoff; E_out[cens] <- 0L
    }
  }
  list(dt = data.table(T = T_out, E = E_out), amb = amb, evt_int = evt_int)
}

# ----------------------------------------------------------------------------
# p_true_i()
# Purpose: the TRUE per-individual confirmation probability for a censored
#          first-positive, derived from the data-generating model. This is the
#          oracle weight for the clone-weighted estimator and works for a first
#          positive at ANY interval (including test 1). See
#          docs/true_p_i_derivation.md for the full derivation.
#
#          Closed form (exponential onset, sensitivity = 1, false-positive
#          rate eps), in terms of the two surrounding test gaps g1, g2:
#            a  = 1 - exp(-lambda*g1)   # P(onset during gap 1 | condition-free)
#            b  = exp(-lambda*g1)       # P(survive gap 1)
#            cc = 1 - exp(-lambda*g2)   # P(onset during gap 2 | condition-free)
#            d  = exp(-lambda*g2)       # P(survive gap 2)
#            p  = (a + eps*b*cc + eps^2*b*d) / (a + eps*b)
# Args:
#   O        - n x K matrix of test (observation) times.
#   interval - integer vector of length n giving each individual's first-
#              positive interval i (e.g. compute_double()$amb); NA where there
#              is no censored first-positive (those entries return NA).
#   lambda   - rate of the exponential onset, T0 ~ Exp(lambda).
#   gap_next - mean gap to the confirming (next) test. The confirming test
#              follows a positive, so this is the after-positive gap mean.
#   eps      - false-positive rate.
# Returns: numeric vector of length n. For an ambiguous individual at interval i
#          it uses the OBSERVED gap g1 = O[, i] - O[, i-1] (with O[, 0] := 0 at
#          the time origin, so interval 1 uses g1 = O[, 1]) and the EXPECTED
#          confirming gap g2 = gap_next (the next test time is unobserved for a
#          censored first-positive, hence the expectation). NA elsewhere.
p_true_i <- function(O, interval, lambda, gap_next, eps) {
  out <- rep(NA_real_, nrow(O))
  idx <- which(!is.na(interval))
  if (length(idx) == 0L) return(out)

  Oaug   <- cbind(0, O)               # col 1 = time origin; col j+1 holds O[, j]
  i      <- interval[idx]
  O_i    <- Oaug[cbind(idx, i + 1L)]  # first-positive test time O[, i]
  O_prev <- Oaug[cbind(idx, i)]       # previous test time O[, i-1] (0 if i == 1)

  g1 <- O_i - O_prev             # observed gap into the first positive
  g2 <- gap_next                 # expected confirming gap (after a positive)
  a  <- 1 - exp(-lambda * g1)    # P(onset during gap 1 | condition-free)
  b  <- exp(-lambda * g1)        # P(survive gap 1 condition-free)
  cc <- 1 - exp(-lambda * g2)    # P(onset during gap 2 | condition-free)
  d  <- exp(-lambda * g2)        # P(survive gap 2 condition-free)
  out[idx] <- (a + eps * b * cc + eps^2 * b * d) / (a + eps * b)
  out
}

# ----------------------------------------------------------------------------
# p_hat_empirical()
# Purpose: a simple data-driven estimate of the (global, scalar) confirmation
#          probability from the test-1 -> test-2 transition. Provided as an
#          alternative to p_true_i() for users who must estimate p from data
#          rather than knowing it. See docs/true_p_i_derivation.md for caveats
#          (single global value, MAR-style exchangeability, horizon selection).
# Args:
#   O      - n x K matrix of test times.
#   RE     - n x K matrix of observed 0/1 results.
#   cutoff - analysis horizon.
# Returns: scalar estimate of p --
#          P(positive at tests 1 AND 2) / P(positive at test 1), restricted to
#          individuals whose second test falls within the horizon.
p_hat_empirical <- function(O, RE, cutoff) {
  mean(RE[, 2] == 1L & RE[, 1] == 1L & O[, 2] < cutoff) /
    mean(RE[, 1] == 1L              & O[, 2] < cutoff)
}

# ----------------------------------------------------------------------------
# build_clone_weighted_data()
# Purpose: expand the double-positive dataset into a clone-weighted one. Each
#          censored first-positive is split into two weighted "clones": an event
#          clone (weight p) and a censored clone (weight 1 - p). A standard
#          weighted KM on the result recovers the omitted events in expectation.
#          Clones censored first-positives at ANY interval (test 1 included),
#          using each individual's own first-positive interval in `amb`.
# Args:
#   O       - n x K matrix of test (observation) times.
#   doub_dt - data.table(T, E) from compute_double() (one row per person).
#   amb     - integer vector (compute_double()$amb): each individual's first-
#             positive interval, or NA if not a censored first-positive.
#   p       - confirmation probability/probabilities supplied by the caller, so
#             anyone can plug in their own weights. Accepts ANY of:
#               * a single scalar (recycled to all ambiguous individuals),
#               * a length-n vector (per individual; subset internally), or
#               * a length-(#ambiguous) vector (already in ambiguous order).
#             Use p_true_i() for the oracle weight, or p_hat_empirical() (or
#             your own estimate) for a data-driven one.
# Returns: data.table(T, E, weights), ready for survfit(..., weights = weights).
build_clone_weighted_data <- function(O, doub_dt, amb, p) {
  # Start from the double-positive data; every individual gets weight 1.
  base <- copy(doub_dt)
  base[, weights := 1.0]

  # Ambiguous individuals = those with a recorded censored-first-positive
  # interval; if none, there is nothing to clone -- return base unchanged.
  is_amb <- !is.na(amb)
  if (!any(is_amb)) return(base)
  n_amb   <- sum(is_amb)
  amb_idx <- which(is_amb)

  # Resolve p to one confirmation probability per ambiguous individual, so the
  # caller may pass a scalar, a full length-n vector, or a pre-subset vector.
  p_amb <-
    if (length(p) == 1L)           rep(p, n_amb)
    else if (length(p) == nrow(O)) p[is_amb]
    else if (length(p) == n_amb)   p
    else stop("`p` must have length 1, nrow(O), or the number of ambiguous individuals.")

  # Each ambiguous individual's event time = their own first-positive test time
  # O[j, amb[j]] (interval-aware, so test-1 cases use O[, 1], test-2 use O[, 2]).
  evt_time <- O[cbind(amb_idx, amb[amb_idx])]

  amb_rows <- base[is_amb]

  # Event clone: the first positive WOULD have been confirmed. Set the event
  # time to that first-positive time, mark E = 1, weight by p.
  cl_event <- copy(amb_rows)
  cl_event[, T       := evt_time]
  cl_event[, E       := 1L]
  cl_event[, weights := p_amb]

  # Censored clone: the first positive would NOT have been confirmed. Keep the
  # original censored (T, E) from doub_dt, weight by the complement 1 - p.
  cl_cens <- copy(amb_rows)
  cl_cens[, weights := 1 - p_amb]

  # Reassemble: the unambiguous rows (weight 1) plus the two clones for each
  # ambiguous individual (weights p and 1 - p, which sum to 1).
  rbind(base[!is_amb], cl_event, cl_cens)
}

# ----------------------------------------------------------------------------
# S_at()
# Purpose: read the survival probability S(t) at a single time t off a fitted
#          KM curve.
# Args:
#   km - a survfit object.
#   t  - the time at which to evaluate S.
# Returns: scalar S(t); extend = TRUE allows t past the last event time, and
#          NA is returned if the curve yields no value.
S_at <- function(km, t) {
  s <- summary(km, times = t, extend = TRUE)$surv
  if (length(s) == 0L) NA_real_ else s
}

# ----------------------------------------------------------------------------
# run_one_rep()
# Purpose: run one full simulation replicate -- generate data, fit the KM
#          estimators, and read each survival curve at the cutoff.
# Args: the same generating parameters as simulate_data(), plus
#   cutoff  - the time at which the curves are compared.
# Returns: a one-row data.table with S(cutoff) for each estimator --
#   True     - error-free single-positive (gold standard).
#   Single   - naive single-positive on observed results.
#   Double   - naive two-consecutive-positives, censored first-positives dropped.
#   Clone    - proposed clone-weighted estimator (oracle weight from p_true_i()).
run_one_rep <- function(n, lambda, n_tests, gap_pos, gap_neg, eps, cutoff) {
  d <- simulate_data(n, lambda, n_tests, gap_pos, gap_neg, eps)

  km_true <- survfit(Surv(T, E) ~ 1, data = compute_single(d$O, d$R,  cutoff))
  km_sing <- survfit(Surv(T, E) ~ 1, data = compute_single(d$O, d$RE, cutoff))

  doub <- compute_double(d$O, d$RE, cutoff)
  km_doub <- survfit(Surv(T, E) ~ 1, data = doub$dt)

  p_vec <- p_true_i(d$O, doub$amb, lambda, gap_pos, eps)   # oracle confirmation prob
  cw_dt <- build_clone_weighted_data(d$O, doub$dt, doub$amb, p_vec)
  km_cw <- survfit(Surv(T, E) ~ 1, data = cw_dt, weights = weights)

  out <- data.table(
    True   = S_at(km_true, cutoff),
    Single = S_at(km_sing, cutoff),
    Double = S_at(km_doub, cutoff),
    Clone  = S_at(km_cw,   cutoff)
  )
  out
}

# -- Style palettes + bias plotter -------------------------------------------
col_pal <- c(True = "black", Single = "blue",
             Double = "red", Clone  = "darkgreen")
sym_pal <- c(Single = 16, Double = 17, Clone = 15)   # circle / triangle / square
lty_pal <- c(True = 1, Single = 2, Double = 2, Clone = 2)
est_label <- c(Single = "Single-positive", Double = "Double-positive",
               Clone  = "Clone-weighted")

# ----------------------------------------------------------------------------
# plot_bias()
# Purpose: scatter+line plot of mean bias vs a swept parameter, one series per
#          estimator. The series drawn are controlled by `ests`.
# Args:
#   dt_long   - long-format data.table with columns Estimator, Bias, and x_var.
#   x_var     - name of the x-axis column ("eps" or "gap_pos").
#   x_label   - x-axis label.
#   file_name - output PDF name (written to out_dir when SAVE_PLOTS).
#   ests      - character vector of estimators to draw (default the main three).
plot_bias <- function(dt_long, x_var, x_label, file_name,
                      ests = c("Single", "Double", "Clone")) {
  n_est    <- length(ests)
  ncol_leg <- if (n_est <= 3) n_est else 2     # 3 series -> one row; 4 -> 2 x 2
  n_row    <- ceiling(n_est / ncol_leg)

  if (SAVE_PLOTS) pdf(file.path(out_dir, file_name),
                      width = 6, height = 5)
  # Enlarge the bottom margin to hold the legend below the plot, then restore.
  op <- par(mar = c(5 + 2 * n_row, 4, 2, 2) + 0.1)
  on.exit(par(op), add = TRUE)

  ylim <- range(dt_long$Bias, na.rm = TRUE)
  ylim <- c(min(0, ylim[1]) - 0.02, max(0, ylim[2]) + 0.02)
  xs   <- dt_long[[x_var]]

  plot(xs, dt_long$Bias, type = "n", ylim = ylim, xaxt = "n",
       xlab = x_label,
       ylab = sprintf("Bias of event rate at cutoff t = %g", cutoff))
  # x-axis ticks exactly at the swept grid values (not R's default pretty ticks).
  xticks <- sort(unique(xs))
  axis(1, at = xticks, labels = formatC(xticks, format = "f", digits = 2))
  abline(h = 0, col = "grey60", lty = 2)

  for (est in ests) {
    sub <- dt_long[Estimator == est]
    sub <- sub[order(sub[[x_var]])]
    points(sub[[x_var]], sub$Bias,
           pch = sym_pal[est], col = col_pal[est], cex = 1.5)
    lines (sub[[x_var]], sub$Bias,
           col = col_pal[est], lwd = 1.5)
  }

  # Legend OUTSIDE the data region: pinned to the bottom-centre of the figure
  # (below the x-axis label), extending up into the enlarged bottom margin.
  legend(x = grconvertX(0.5, from = "ndc", to = "user"),
         y = grconvertY(0.02, from = "ndc", to = "user"),
         xjust = 0.5, yjust = 0, xpd = NA,
         legend = est_label[ests],
         pch    = sym_pal[ests],
         col    = col_pal[ests],
         lty = 1, lwd = 1.5, bty = "n", ncol = ncol_leg, cex = 0.85)

  if (SAVE_PLOTS) dev.off()
}
```

### 3.2 `Simulation.R` — main figures

```r
# ============================================================================
# Simulation -- main figures
#
# Sources scripts/sim_core.R (shared functions + parameters) and produces the
# three main poster figures:
#   * sim_main_km.pdf      -- KM curves at baseline (True/Single/Double/Clone)
#   * sim_bias_vs_eps.pdf  -- bias vs false-positive rate (epsilon)
#   * sim_bias_vs_gap.pdf  -- bias vs after-positive gap mean
#
# Sweep design (one parameter varied at a time, the other held at baseline):
#   Sweep A: epsilon in {0.02, 0.05, 0.10, 0.15}, after-positive gap = 0.4
#   Sweep B: after-positive gap in {0.25, 0.40, 0.55, 0.75}, epsilon = 0.05
#            (after-negative gap held at 0.75)
# ============================================================================

source(if (file.exists("sim_core.R")) "sim_core.R" else "scripts/sim_core.R")

# -- Run sweeps --------------------------------------------------------------

cat(sprintf("Running epsilon sweep (%d settings x %d reps)...\n",
            length(eps_grid), n_reps))
sweep_eps <- rbindlist(lapply(eps_grid, function(eps) {
  reps <- rbindlist(lapply(seq_len(n_reps), function(r)
    run_one_rep(n, lambda1, n_tests,
                gap_pos_baseline, gap_neg_baseline, eps, cutoff)))
  reps[, eps := eps]
  reps
}))

cat(sprintf("Running gap sweep (%d settings x %d reps)...\n",
            length(gap_grid), n_reps))
sweep_gap <- rbindlist(lapply(gap_grid, function(gp) {
  reps <- rbindlist(lapply(seq_len(n_reps), function(r)
    run_one_rep(n, lambda1, n_tests, gp, gap_neg_baseline, eps_baseline, cutoff)))
  reps[, gap_pos := gp]
  reps
}))

# -- Aggregate to mean bias on event rate at cutoff --------------------------
# bias_estimator = (1 - S_est) - (1 - S_true) = S_true - S_est
# (signed: positive = estimator overstates event rate; negative = understates)

bias_eps <- sweep_eps[, .(
  Single = mean(True - Single, na.rm = TRUE),
  Double = mean(True - Double, na.rm = TRUE),
  Clone  = mean(True - Clone,  na.rm = TRUE)
), by = eps]

bias_gap <- sweep_gap[, .(
  Single = mean(True - Single, na.rm = TRUE),
  Double = mean(True - Double, na.rm = TRUE),
  Clone  = mean(True - Clone,  na.rm = TRUE)
), by = gap_pos]

cat("\nBias vs epsilon:\n");                print(bias_eps)
cat("\nBias vs after-positive gap:\n");     print(bias_gap)

bias_eps_long <- melt(bias_eps, id.vars = "eps",
                      variable.name = "Estimator", value.name = "Bias")
bias_gap_long <- melt(bias_gap, id.vars = "gap_pos",
                      variable.name = "Estimator", value.name = "Bias")

# -- Plot 1: Main KM curves at baseline --------------------------------------
cat(sprintf("\nGenerating main KM plot at baseline (n = %d)...\n", n))
d0 <- simulate_data(n, lambda1, n_tests,
                    gap_pos_baseline, gap_neg_baseline, eps_baseline)

km_true <- survfit(Surv(T, E) ~ 1, data = compute_single(d0$O, d0$R,  cutoff))
km_sing <- survfit(Surv(T, E) ~ 1, data = compute_single(d0$O, d0$RE, cutoff))
doub0   <- compute_double(d0$O, d0$RE, cutoff)
km_doub <- survfit(Surv(T, E) ~ 1, data = doub0$dt)
p_vec0  <- p_true_i(d0$O, doub0$amb, lambda1, gap_pos_baseline, eps_baseline)
cw_dt0  <- build_clone_weighted_data(d0$O, doub0$dt, doub0$amb, p_vec0)
km_cw   <- survfit(Surv(T, E) ~ 1, data = cw_dt0, weights = weights)

km_ests <- c("True", "Single", "Double", "Clone")
if (SAVE_PLOTS) pdf(file.path(out_dir, "sim_main_km.pdf"),
                    width = 6, height = 4.5)
plot(km_true, conf.int = FALSE, xlim = c(0, cutoff), ylim = c(0, 1),
     xlab = "Time", ylab = "Probability of remaining undiagnosed",
     col = col_pal["True"], lty = lty_pal["True"], lwd = 2)
lines(km_sing, conf.int = FALSE,
      col = col_pal["Single"], lty = lty_pal["Single"], lwd = 2)
lines(km_doub, conf.int = FALSE,
      col = col_pal["Double"], lty = lty_pal["Double"], lwd = 2)
lines(km_cw,   conf.int = FALSE,
      col = col_pal["Clone"],  lty = lty_pal["Clone"],  lwd = 2)
legend("bottomleft",
       legend = c("True", "Single-positive", "Double-positive", "Clone-weighted"),
       col    = col_pal[km_ests],
       lty    = lty_pal[km_ests],
       lwd = 2, bty = "n")
if (SAVE_PLOTS) dev.off()

# -- Plot 2 / Plot 3: Bias scaling plots -------------------------------------
cat("Generating bias-vs-epsilon plot...\n")
plot_bias(bias_eps_long, "eps",
          "False-positive rate (epsilon)",
          "sim_bias_vs_eps.pdf")

cat("Generating bias-vs-gap plot...\n")
plot_bias(bias_gap_long, "gap_pos",
          "Mean test gap after a positive result",
          "sim_bias_vs_gap.pdf")

cat("\nDone.\n")
if (SAVE_PLOTS) cat("Figures saved to:", normalizePath(out_dir), "\n")
```

---

## 4. Inference for the proposed estimator

> **Work in progress.** The approach presented here outlines a possible route to
> inference for the proposed estimator. It is not intended as a rigorous proof.
> A more rigorous derivation of the inference is in development.

The poster reports point estimates and Monte-Carlo bias. This section sets out
the notation, derives the estimator from the standard Kaplan–Meier formula, and
briefly describes where inference stands.

### 4.1 Setting and notation

The notation is written for a general **trigger and confirmation** structure, so
that the argument does not depend on the diagnostic-testing example. An event of
interest is recorded only when a provisional first event is followed by a
confirming one. Notation is local to this section.

| Symbol | Meaning | In the running example |
| :-- | :-- | :-- |
| $X_0$ | time of the **event of interest**, a trigger that is confirmed | time of the first of two consecutive positive tests |
| $X_1$ | time of the **trigger event**, the provisional event that opens a confirmation | time of the first positive test |
| $X_2$ | time of the **confirmation event** that follows the trigger, that is, the time at which $Y$ is measured | time of the confirmatory test |
| $Y$ | binary result of the confirmation event: $Y = 1$ confirmed, $Y = 0$ not confirmed | result of the confirmatory test, positive or negative |
| $C$ | censoring time | administrative end of follow-up at the analysis horizon |

Individuals are indexed by $i = 1, \dots, n$, and the index is suppressed
wherever only one individual is in view. The event of interest is the trigger
that the confirmation upholds, so $X_{0,i}$ is determined by $X_{1,i}$ and $Y_i$:

$$
X_{0,i} = X_{1,i} \ \text{ if } Y_i = 1, \qquad
X_{0,i} = \infty \ \text{ if } Y_i = 0 ,
$$

where $X_{0,i} = \infty$ records that no event of interest occurs.

A confirmation occurs after a trigger that has already taken place, so by
definition

$$
X_{1,i} < X_{2,i} .
$$

An individual may pass through several trigger events over follow-up. A refuted
trigger does not end observation, so a later trigger can still arise and be
confirmed, and the event of interest is the first trigger that is confirmed. For
simplicity the notation does not carry this complexity, and what follows is
written for a single trigger and its confirmation.

Let $R_i$ indicate the ambiguous case, in which the trigger is observed but the
confirmation is not, in other words, we cannot verify whether
$X_{0,i} = X_{1,i}$:

$$
R_i = \mathbf{1}\lbrace X_{1,i} \le C_i < X_{2,i}\rbrace .
$$

Take first the pair that an ordinary survival analysis would use:

$$
Z_{0,i} = \min(X_{0,i}, C_i), \qquad
\Delta_{0,i} = \mathbf{1}\lbrace X_{0,i} \le C_i\rbrace .
$$


In our scenario this pair cannot be observed for every individual, because
$X_{0,i}$ itself may be unknown. The subscript $0$ is deliberate: it marks the
quantities we would like to have, as against the data actually observed.

$Z_{0,i}$ and $\Delta_{0,i}$ are what the analysis has to recover from the data.
Once they are in hand, they can be entered into the ordinary Kaplan–Meier
formula, which estimates the estimand of interest,

$$
S(t) = P(X_0 > t) .
$$

and can be rewritten as

$$
S(t) = \prod_{k \le t}
\left\lbrace 1 - \frac{P[Z_0 = k, \ \Delta_0 = 1]}
{P[Z_0 \ge k]} \right\rbrace . \qquad (1)
$$

Everything in §4.2 is aimed at the two probabilities in equation (1).

The pair $(Z_{0,i}, \Delta_{0,i})$ is available only when $R_i = 0$. When
$R_i = 1$ the value of $Y_i$ is never seen, so $X_{0,i}$ is unknown and the pair
is one of two possibilities: $(X_{1,i}, 1)$ if $Y_i = 1$, and $(C_i, 0)$ if
$Y_i = 0$.

Even though $X_{0,i}$ is not observed in the ambiguous case, we know how it is
determined by the unobserved $Y_i$. The pair can therefore be written out across
the two cases:

$$
Z_{0,i} = (1 - R_i)\min(X_{0,i}, C_i) +
R_i \mathbf{1}\lbrace Y_i = 1\rbrace X_{1,i} +
R_i \mathbf{1}\lbrace Y_i = 0\rbrace C_i , \qquad (2)
$$

$$
\Delta_{0,i} = (1 - R_i)\mathbf{1}\lbrace X_{0,i} \le C_i\rbrace +
R_i \mathbf{1}\lbrace Y_i = 1\rbrace . \qquad (3)
$$

Everything on the right of these two expressions is observed except $Y_i$.

For individual $i$ observed data are

$$
O_i = (R_i, (1 - R_i) Z_{0,i}, (1 - R_i) \Delta_{0,i},
 R_i X_{1,i}, R_i C_i) .
$$

For a non-ambiguous individual, $R_i = 0$, this reduces to the usual survival
pair $(Z_{0,i}, \Delta_{0,i})$. For an ambiguous one, $R_i = 1$, it holds the
trigger time $X_{1,i}$ and the censoring time $C_i$, with $X_{0,i}$
missing.

The case $R_i = 1$ is the one this work addresses. It is an *ambiguity* rather
than sampling uncertainty: nothing about the individual is noisy, and the single
datum that would settle the event status falls outside the observation window.
Clone weighting resolves it by splitting the individual into an event clone of
weight $p_i$ and a censored clone of weight $1 - p_i$. The next subsection makes
that precise and defines $p_i$.

### 4.2 Estimator derivation

Equation (1) needs two population quantities at each event time $k$. Take the
denominator first. The Kaplan–Meier estimator counts the individuals still at
risk, $\sum_i \mathbf{1}\lbrace Z_{0,i} \ge k\rbrace$, whose population counterpart is
$P(Z_0 \ge k)$. Substituting (2) and splitting on the three cases gives

$$
P(Z_0 \ge k) = P[R = 0, \ \min(X_0, C) \ge k] +
P[R = 1, \ Y = 1, \ X_1 \ge k] +
P[R = 1, \ Y = 0, \ C \ge k] .
$$

The three events are disjoint and cover every individual. The first term is computable from the observed data. The other two are
not, because they condition on $Y$, which is never seen when $R = 1$.

Both ambiguous terms can be written exactly by
conditioning on what is observed. For an individual with $R_i = 1$, let

$$
p_i = P(Y_i = 1 \mid O_i)
$$

be the probability that their trigger would have been confirmed, given
everything observed about them. Iterated expectation then gives

$$
P[R = 1, \ Y = 1, \ X_1 \ge k]
 = E[R \mathbf{1}\lbrace X_1 \ge k\rbrace p_i] ,
$$

$$
P[R = 1, \ Y = 0, \ C \ge k]
 = E[R \mathbf{1}\lbrace C \ge k\rbrace (1 - p_i)] .
$$

It simply means walking through the ambiguous individuals still at risk at
time $k$ and adding up their individual weights.

Substituting both back into the at-risk probability gives

$$
P(Z_0 \ge k) = P[R = 0, \ \min(X_0, C) \ge k] +
E[R \mathbf{1}\lbrace X_1 \ge k\rbrace p_i] +
E[R \mathbf{1}\lbrace C \ge k\rbrace (1 - p_i)] , \qquad (4)
$$

which is what the two clones compute: an event clone placed at $X_{1,i}$ and
carrying weight $p_i$, and a censored clone placed at $C_i$ and carrying weight
$1 - p_i$.

The numerator follows the same way, from (3). Splitting
$\lbrace Z_0 = k, \ \Delta_0 = 1\rbrace$ on the same three cases, an ambiguous individual
with $Y = 0$ has $\Delta_0 = 0$ and cannot contribute, so only two terms remain:

$$
P[Z_0 = k, \ \Delta_0 = 1]
 = P[R = 0, \ X_0 = k \le C] +
E[R \mathbf{1}\lbrace X_1 = k\rbrace p_i] . \qquad (5)
$$

So $p_i$ enters both the numerator and the denominator, while $1 - p_i$ enters
the denominator only. That asymmetry is the mechanism: the censored clone extends
time at risk without producing an event.

The cumulative version, with $Z_0 \le k$ in place of $Z_0 = k$, is derived in the
same way. It is the form to use when event times are continuous rather than
confined to the inspection grid.

In both (4) and (5) the only ingredient that is not identifiable from the
observed data is $p_i$. Everything else is a function of $O_i$. So if $p_i$ is
known, or can be estimated consistently, substituting it gives consistent
estimators of $P[Z_0 \ge k]$ and $P[Z_0 = k, \Delta_0 = 1]$, and through (1) a
consistent estimator of $S(t)$.

Finally, given $O_i$, the two quantities can be estimated directly. Writing
$\hat p_i$ for the known or estimated weight, their sample counterparts are

$$
\widehat{P}[Z_0 \ge k]
 = \frac{1}{n}\sum_{i=1}^{n} \Big[(1 - R_i)\mathbf{1}\lbrace Z_{0,i} \ge k\rbrace +
R_i \hat p_i \mathbf{1}\lbrace X_{1,i} \ge k\rbrace +
R_i (1 - \hat p_i) \mathbf{1}\lbrace C_i \ge k\rbrace
\Big] , \qquad (6)
$$

$$
\widehat{P}[Z_0 = k, \ \Delta_0 = 1]
 = \frac{1}{n}\sum_{i=1}^{n} \Big[(1 - R_i)\mathbf{1}\lbrace Z_{0,i} = k, \ \Delta_{0,i} = 1\rbrace +
R_i \hat p_i \mathbf{1}\lbrace X_{1,i} = k\rbrace
\Big] . \qquad (7)
$$

Every term is computable from $O_i$. The first summand is the pair that a
non-ambiguous individual supplies directly, and the other two are the trigger and
censoring times that an ambiguous individual supplies, carrying weights
$\hat p_i$ and $1 - \hat p_i$. Substituting (6) and (7) into (1) gives the
clone-weighted estimator of $S(t)$.

This is exactly the cloning construction. Read (6) and (7) as a weighted
Kaplan–Meier on a dataset holding one row per non-ambiguous individual, at
$(Z_{0,i}, \Delta_{0,i})$ with weight 1, and two rows per ambiguous individual,
an event clone at $(X_{1,i}, 1)$ with weight $\hat p_i$ and a censored clone at
$(C_i, 0)$ with weight $1 - \hat p_i$. That is the dataset
`build_clone_weighted_data()` assembles in §3, and passing it to a weighted
Kaplan–Meier routine evaluates (6) and (7) term by term.

The two clones do not double count. Their weights sum to one, so up to the
trigger time both are at risk and the individual contributes exactly one person
to the denominator, the same as anybody else. The split takes effect only after
$X_{1,i}$: a fraction $\hat p_i$ leaves the risk set as an event, and the
remaining $1 - \hat p_i$ stays at risk until $C_i$. Two rows are needed because an ambiguous
individual has two possible exit times as well as two possible event statuses,
and a single row can hold only one exit time.

The same expression also contains the naive double-positive estimator of §2.4,
the comparator used in the simulation, as the case $\hat p_i = 0$. Every censored
first-positive is then left as a non-event censored at $C_i$. The formula shows
directly why that understates the event rate. The numerator (7) loses the term
$R_i \hat p_i \mathbf{1}\lbrace X_{1,i} = k\rbrace$, so those events are never counted,
while the denominator (6) keeps the same individuals at risk with full weight all
the way to $C_i$. Fewer events over an undiminished risk set gives a smaller
hazard at every $k$, so the survival curve sits too high. That is the negative
bias of the double-positive curve in Figure 2 of the poster.

Where $p_i$ comes from is the subject of §5.

### 4.3 Inference

Inference for this estimator is a topic for future work. What follows records the
routes being considered rather than settled results.

The unit of resampling should be the individual. A cluster bootstrap that
resamples individuals and rebuilds their clones within each replicate accounts
for the dependence between an individual's two clones.

How the weight is treated then depends on where it came from, following the three
routes of §5.

- **$p_i$ known.** Treat it as a constant and take bootstrap intervals with the
  weight plugged in. The only randomness is the sampling of individuals.
- **$p_i$ taken from an external source.** Draw that source's parameters once per
  bootstrap replicate, according to their reported variation, and recompute every
  $p_i$ from the draw. The interval then carries the uncertainty in the external
  estimate as well as the sampling variability of the study data.
- **$p_i$ estimated from the study data.** Embed the estimation of $p_i$ inside
  the bootstrap, refitting it on every resampled dataset, so that the variability
  of the weight is propagated rather than held fixed.

---

## 5. Estimating *pᵢ*

Whether *pᵢ* is known or has to be estimated depends on the setting. In some
applications the confirmation step is specified, and *pᵢ* follows from that
specification.

Equipment fault monitoring is one such case. A fault is often logged only when an
alarm repeats on the next inspection cycle, which is similar to the rule of two
consecutive test results used here. The alarm rate, whether genuine or false, is
fixed by the sensor design and published by the manufacturer. Both parts of the
confirmation probability are known before any data are collected.

On the other hand, an external data source can provide useful information to
derive *pᵢ*. A validation sample may report the operating characteristics of the
test directly: sensitivity, specificity, and the confirmation rate itself.
Published prevalence from disease registries can provide true disease status by
time and tell how often a first positive is genuine rather than spurious. These
sources often approximate *pᵢ* closely enough that no model is needed. As a last
resort, *pᵢ* can be estimated from a model fitted to the observed test results in
the study at hand.

In the simulation of §2, *pᵢ* is computed exactly from the closed form derived in
[`docs/true_p_i_derivation.md`](docs/true_p_i_derivation.md), using the
false-positive rate ε and the onset rate λ that generated the data; both are
known in the simulation. This simplification is deliberate to focus on how
clone-based weighting works, so the weight is held exact and the separate problem
of estimating *pᵢ* is put aside. Any residual bias in the clone-weighted curve of
Figure 2 is therefore a property of the event definition, discussed in §6, and
not of a misspecified *pᵢ*.

---

## 6. Residual false-positive bias and its correction

This is the discussion referenced in Figure 2 of the poster, where the
clone-weighted curve does not sit exactly at zero and its bias grows with the
false-positive rate ε.

**Why the residual bias is there.** Clone weighting resolves an *ambiguity* about
censored first-positive cases: a first positive whose following confirming test
is never observed is neither a confirmed event nor a confirmed non-event. What
clone weighting does not do is resolve the defects inherent in the event
*definition* itself (the confirmation rule). Under this rule, a pair of
consecutive positives is treated as a confirmed diagnosis even when both
positives are false, which happens with probability of order ε² per pair of
tests. The double-positive definition is defective in its own way, and any
estimator targeting this definition, including the clone-weighted one, has this
bias. This is why, in Figure 2 of the poster, the residual bias grows with ε
while the bias from censored first-positives is removed.

The same mechanism explains a feature of the naive double-positive estimator that
might otherwise be mistaken for an improvement, as its bias decreases with ε. Two
errors are at work in opposite directions. Discarding censored first-positives
underestimates the event rate, spurious confirmations overestimate it, and only
the second grows with ε. At larger ε the two nearly cancel, so the
double-positive estimator lands closer to the truth, not because it is unbiased
but because its two sources of error offset each other. Clone weighting removes
the first error but leaves the second one exposed.

**How it can be handled.** Specificity is routinely reported for diagnostic
assays, so ε = 1 − specificity is often known or can be well approximated. When
it is, the expected number of spurious confirmations is a function of ε and the
observed test results, and it can be removed from the event counts so that they
reflect only true positives. This belongs to the broader problem of outcome
misclassification in time-to-event analysis, which has a rich literature of its
own. When ε is not known externally, a validated subsample is the usual route.

A detailed discussion of corrections for mismeasured outcomes, and of how they
interact with the clone weights, is a separate line of work and is deliberately
left out of this supplement.

**Further reading**

- Meier AS, Richardson BA, Hughes JP. Discrete proportional hazards models for
  mismeasured outcomes. *Biometrics.* 2003;59(4):947–954.
- Balasubramanian R, Lagakos SW. Estimation of a failure time distribution based
  on imperfect diagnostic tests. *Biometrika.* 2003;90(1):171–182.

---

**Use of AI tools.** Claude (Anthropic) was used in preparing this repository.
The conceptual work is the author's throughout, from posing the research question
to proposing the method that answers it. The tool assisted by expanding and
critiquing those ideas, editing prose, checking derivations, and writing and
documenting the simulation code. All modelling choices and final wording are the
author's, who reviewed and is responsible for everything here.
