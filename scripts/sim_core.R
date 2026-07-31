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
