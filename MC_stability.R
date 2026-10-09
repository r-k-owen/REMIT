## ============================================================================
## MC_stability.R
## Assess Monte Carlo stability of the patient-level simulation (PLS) 
## (B = 100).
##
## Inputs (written by Patient_level_simulation.R in the working directory):
##   pooled_overall_incl_eligible.csv  - one row per iteration: overall and
##                                       treatment-eligible life expectancy
##                                       differences (SGLT-2i minus SOC)
##   per_iteration_differences.csv     - one row per iteration x transition:
##                                       paired differences in transition
##                                       probability, age at transition and
##                                       life expectancy
##
## Outputs:
##   MC_stability_summary.csv   - estimate, 95% interval, Monte Carlo SEs and
##                                the number of iterations needed to meet the
##                                chosen tolerance, for every outcome
##   MC_trace_key_outcomes.png  - running mean and running 2.5th/97.5th
##                                percentiles against iteration number
##   MC_trace_all_LE.png        - as above for every life expectancy difference
##
## ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(2026)        # for the bootstrap only
n_boot <- 2000        # bootstrap resamples of the B iteration-level results

## Tolerances for the MCSE, on the scale results are reported (2 d.p. for
## years, 3 d.p. for probabilities). 
tol_years <- 0.05     # life expectancy and age at transition, years
tol_prob  <- 0.005    # transition probabilities

## ---- Read inputs ----
overall <- read.csv("pooled_overall_incl_eligible.csv")
diffs   <- read.csv("per_iteration_differences.csv")

## ---- Relabel K/D/H state codes as CKD/T2DM/HF ----
relabel_state <- function(x) {
  map <- c(K = "CKD", D = "T2DM", H = "HF")
  vapply(strsplit(x, " ", fixed = TRUE), function(words) {
    paste(vapply(words, function(w) {
      if (grepl("^[KDH]{1,3}$", w)) paste(map[strsplit(w, "")[[1]]], collapse = "-") else w
    }, character(1)), collapse = " ")
  }, character(1))
}
diffs$label <- relabel_state(diffs$label)

## ---- Long format: one row per outcome x iteration ----
long <- bind_rows(
  overall %>% transmute(outcome = "Life expectancy: treatment-eligible population",
                        measure = "Life expectancy", iter, value = diff_life_eligible),
  overall %>% transmute(outcome = "Life expectancy: whole simulated population",
                        measure = "Life expectancy", iter, value = diff_life),
  diffs %>% transmute(outcome = paste0(label, ": life expectancy"),
                      measure = "Life expectancy", iter, value = life_diff),
  diffs %>% transmute(outcome = paste0(label, ": transition probability"),
                      measure = "Transition probability", iter, value = prob_diff),
  diffs %>% transmute(outcome = paste0(label, ": age at transition"),
                      measure = "Age at transition", iter, value = age_diff)
) %>%
  filter(is.finite(value)) %>%
  arrange(outcome, iter)

## Transitions from triple-condition states to death have probability 1.
long <- long %>%
  group_by(outcome) %>%
  filter(!(measure == "Transition probability" & sd(value) == 0)) %>%
  ungroup()

## ---- Monte Carlo standard errors and required iterations ----
mc_summary <- long %>%
  group_by(outcome, measure) %>%
  group_modify(function(d, key) {
    v <- d$value
    B <- length(v)
    boots <- replicate(n_boot, {
      s <- sample(v, B, replace = TRUE)
      quantile(s, c(0.025, 0.975), names = FALSE)
    })
    tol <- if (key$measure == "Transition probability") tol_prob else tol_years
    mcse_mean <- sd(v) / sqrt(B)
    mcse_lo   <- sd(boots[1, ])
    mcse_hi   <- sd(boots[2, ])
    half      <- v[seq_len(floor(B / 2))]
    tibble(
      B          = B,
      estimate   = mean(v),
      ci_lo      = quantile(v, 0.025, names = FALSE),
      ci_hi      = quantile(v, 0.975, names = FALSE),
      mcse_mean  = mcse_mean,
      mcse_ci_lo = mcse_lo,
      mcse_ci_hi = mcse_hi,
      tolerance  = tol,
      B_req_mean = ceiling(B * (mcse_mean / tol)^2),
      B_req_ci   = ceiling(B * (max(mcse_lo, mcse_hi) / tol)^2),
      meets_tolerance = mcse_mean <= tol & max(mcse_lo, mcse_hi) <= tol,
      estimate_first_half = mean(half),
      change_half_to_full = mean(v) - mean(half)
    )
  }) %>%
  ungroup() %>%
  arrange(measure, desc(pmax(B_req_mean, B_req_ci)))

write.csv(mc_summary, "MC_stability_summary.csv", row.names = FALSE)

## ---- Console report ----
cat("\nMonte Carlo stability summary (B =", max(mc_summary$B), "iterations)\n")
cat(sprintf("Tolerance: %.3f years for life expectancy / age; %.3f for probabilities\n",
            tol_years, tol_prob))
cat(sprintf("Outcomes meeting tolerance for mean and both interval limits: %d of %d\n",
            sum(mc_summary$meets_tolerance), nrow(mc_summary)))
key <- mc_summary %>% filter(grepl("treatment-eligible", outcome))
cat(sprintf(paste0("Treatment-eligible life expectancy difference: %.2f (95%% CI %.2f, %.2f); ",
                   "MCSE mean %.3f, lower limit %.3f, upper limit %.3f\n"),
            key$estimate, key$ci_lo, key$ci_hi, key$mcse_mean, key$mcse_ci_lo, key$mcse_ci_hi))
worst <- mc_summary %>% filter(!meets_tolerance)
if (nrow(worst) > 0) {
  cat("\nOutcomes NOT meeting tolerance (iterations required shown):\n")
  print(as.data.frame(worst %>% select(outcome, mcse_mean, mcse_ci_lo, mcse_ci_hi, B_req_mean, B_req_ci)),
        digits = 3, row.names = FALSE)
}

## ---- Running statistics for trace plots ----
running <- long %>%
  group_by(outcome, measure) %>%
  mutate(
    k        = row_number(),
    run_mean = cummean(value),
    run_lo   = vapply(k, function(j) quantile(value[1:j], 0.025, names = FALSE), numeric(1)),
    run_hi   = vapply(k, function(j) quantile(value[1:j], 0.975, names = FALSE), numeric(1))
  ) %>%
  ungroup()

trace_plot <- function(dat, title) {
  finals <- dat %>% group_by(outcome) %>% filter(k == max(k))
  ggplot(dat, aes(x = k)) +
    geom_ribbon(aes(ymin = run_lo, ymax = run_hi), fill = "grey85") +
    geom_line(aes(y = run_lo), colour = "grey40", linewidth = 0.3) +
    geom_line(aes(y = run_hi), colour = "grey40", linewidth = 0.3) +
    geom_line(aes(y = run_mean), colour = "black", linewidth = 0.6) +
    geom_hline(data = finals, aes(yintercept = run_mean), linetype = "dashed", colour = "grey30") +
    geom_hline(yintercept = 0, colour = "red", linewidth = 0.3) +
    facet_wrap(~ outcome, scales = "free_y", labeller = label_wrap_gen(40)) +
    labs(title = title,
         x = "Number of Monte Carlo iterations",
         y = "Difference (SGLT-2 inhibitors minus SOC)",
         caption = "Black: running mean; grey band: running 2.5th-97.5th percentiles; dashed: final estimate") +
    theme_bw(base_size = 9)
}

## Key outcomes plot
key_outcomes <- c(
  "Life expectancy: treatment-eligible population",
  "T2DM -> Death: transition probability",
  "CKD -> CKD-T2DM: age at transition",
  "HF -> HF-T2DM: age at transition",
  "CKD-HF -> CKD-HF-T2DM: age at transition",
  "HF-CKD -> HF-CKD-T2DM: age at transition"
)
p_key <- trace_plot(running %>% filter(outcome %in% key_outcomes),
                    "Monte Carlo stability: key outcomes")
ggsave("MC_trace_key_outcomes.png", p_key, width = 10, height = 6, dpi = 300)

p_le <- trace_plot(running %>% filter(measure == "Life expectancy"),
                   "Monte Carlo stability: life expectancy differences")
ggsave("MC_trace_all_LE.png", p_le, width = 14, height = 12, dpi = 300)

cat("\nWritten: MC_stability_summary.csv, MC_trace_key_outcomes.png, MC_trace_all_LE.png\n")
