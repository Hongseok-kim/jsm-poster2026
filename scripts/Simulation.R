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
