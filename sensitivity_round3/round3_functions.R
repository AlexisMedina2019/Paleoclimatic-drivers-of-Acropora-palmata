## Round-3 sensitivity analyses (Section 4.7 of regional_analysis_Final.qmd).
## Sourced inside the qmd: uses the objects of the session (recs, age_draws, build_series,
## cli_tab, formulas_mg, formulas_cb, drop_hydro, count_smooth_terms, cap_smooth_k,
## safe_k_for, term_pvalues, panel_names, best_*_cb).

## ---- crash-safe fit -------------------------------------------------------------------
## On near-binary responses with collinear regimes (e.g. M7 in Punta Maroma) mgcv can abort R
## itself (memory corruption), which tryCatch cannot intercept. On Unix each fit with a
## categorical regime runs in a forked child (killed after 60 s); a crashed fit counts as
## not fitted. On Windows (no fork) fits run directly.
fit_summary <- function(f, d, combined) {
  ns <- count_smooth_terms(f); n <- nrow(d)
  if (n < 8 || n <= ns * 3) return(NULL)
  ff <- if (combined) cap_smooth_k(f, safe_k_for(n, ns)) else f
  m <- tryCatch(gam(ff, data = d, family = betar(link = "logit"), method = "REML", select = combined),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  list(aic = AIC(m), p = tryCatch(term_pvalues(m), error = function(e) numeric(0)))
}
safe_fit <- function(f, d, combined) {
  has_factor <- any(c("Thermal", "Hurricane", "Hydro") %in% all.vars(f))
  if (!has_factor || .Platform$OS.type != "unix") return(fit_summary(f, d, combined))
  j <- parallel::mcparallel(fit_summary(f, d, combined), silent = TRUE)
  r <- parallel::mccollect(j, wait = FALSE, timeout = 60)
  if (is.null(r)) { tools::pskill(j$pid, tools::SIGKILL); parallel::mccollect(j, wait = FALSE); return(NULL) }
  r <- r[[1]]; if (is.null(r) || inherits(r, "try-error")) NULL else r
}
pget <- function(fits, mod, pat) {
  if (is.null(fits[[mod]])) return(NA_real_)
  v <- fits[[mod]]$p[grepl(pat, names(fits[[mod]]$p))]; if (length(v)) min(v) else NA_real_
}

## ---- climate table for a bin origin shifted by o years (bins in years before 1950 + o) ----
## TC event layers (Schmitt et al. 2025, table S4a) re-binned exactly; youngest bin counted from
## the youngest dated colony (1956 CE, -6 yr BP). Temperature interpolated at bin centres;
## regimes from the BP bin that holds the bin centre.
ev_bp  <- read.csv(here("sensitivity_round3", "schmitt2025_tableS4a_event_layers.csv"))$varve_years_before_2022 - 72
cli_bp <- cli_tab                  # BP climate table of the main analysis (Section 4.5.2)
age_ref <- cli_tab$age
make_cli <- function(o) {
  lo <- 100 * (0:59) - o; hi <- lo + 100; ctr <- (lo + hi) / 2
  ib <- pmin(60, pmax(1, floor(ctr / 100) + 1))
  Tn <- approx(seq(50, 5950, 100), cli_bp$GMHT_CPS0_30N, ctr, rule = 2)$y
  hy <- as.character(cli_bp$Hydro_YUC)[ib]
  tibble::tibble(age = age_ref, GMHT_CPS0_30N = Tn, dT_raw = c(0, diff(Tn)),
    Hur100_fqcy = sapply(1:60, function(k) sum(ev_bp >= ifelse(k == 1, -6, lo[k]) & ev_bp < hi[k])),
    Hur100_av   = sapply(1:60, function(k) sum(ev_bp >= lo[k] - 50 & ev_bp < hi[k] + 50) / 2),
    Thermal = as.character(cli_bp$Thermal)[ib], Hydro_YUC = ifelse(is.na(hy), "no data", hy),
    Hydro_FL = as.character(cli_bp$Hydro_FL)[ib], Hur_YUC = as.character(cli_bp$Hur_YUC)[ib],
    Hur_FL = as.character(cli_bp$Hur_FL)[ib],
    Thermal_lag = c(as.character(cli_bp$Thermal)[ib][-1], tail(as.character(cli_bp$Thermal)[ib], 1)))
}
series_at <- function(o) {          # builds the five series with the shifted bins
  old <- cli_tab; on.exit(cli_tab <<- old)
  cli_tab <<- make_cli(o)
  build_series(recs %>% dplyr::select(spot, source, t_ka) %>% dplyr::mutate(t_ka = t_ka + o / 1000))
}

## ---- phase ensemble ----------------------------------------------------------------------
run_phase_ensemble <- function(offsets = seq(0, 90, 10)) {
  dplyr::bind_rows(lapply(offsets, function(o) {
    ser <- series_at(o)
    dplyr::bind_rows(lapply(panel_names, function(p) {
      d  <- ser %>% dplyr::filter(panel == p, !is.na(counts_norm_sv))
      mg <- Filter(Negate(is.null), lapply(drop_hydro(formulas_mg, d), safe_fit, d = d, combined = FALSE))
      cb <- Filter(Negate(is.null), lapply(drop_hydro(formulas_cb, d), safe_fit, d = d, combined = TRUE))
      a_mg <- sapply(mg, `[[`, "aic"); a_cb <- sapply(cb, `[[`, "aic"); best_cb <- names(which.min(a_cb))
      data.frame(o = o, panel = p, n = nrow(d), best_marg = names(which.min(a_mg)),
        p_hf_marg = pget(mg, "M_hf", "Hur100"), p_hav_marg = pget(mg, "M_hav", "Hur100"),
        best_comb = best_cb, dAIC_comb_vs_M1 = unname(a_cb["M1"] - min(a_cb)),
        p_hur_best = pget(cb, best_cb, "Hur100"), p_therm_best = pget(cb, best_cb, "^Thermal$"),
        p_hurreg_best = pget(cb, best_cb, "^Hurricane$"))
    }))
  }))
}
run_phase_regime_contrasts <- function(offsets = seq(0, 90, 10), panels = c("Punta Maroma", "Belize", "Regional")) {
  dplyr::bind_rows(lapply(offsets, function(o) {
    ser <- series_at(o)
    dplyr::bind_rows(lapply(panels, function(p) {
      d <- ser %>% dplyr::filter(panel == p, !is.na(counts_norm_sv))
      d$Hurricane <- relevel(droplevels(factor(d$Hurricane)), ref = "moderate")
      m <- tryCatch(gam(counts_norm_sv ~ s(age) + Hurricane, data = d, family = betar(link = "logit"), method = "REML"),
                    error = function(e) NULL)
      if (is.null(m)) return(NULL)
      pt <- summary(m)$p.table; pt <- pt[grepl("^Hurricane", rownames(pt)), , drop = FALSE]
      data.frame(o = o, panel = p, level = sub("Hurricane", "", rownames(pt)), est = pt[, 1], p = pt[, 4])
    }))
  }))
}

## ---- LIA contrasts without the two modern centuries ------------------------------------
run_lia_tests <- function(offsets = seq(0, 90, 10), panels = c("Punta Maroma", "Belize", "Regional")) {
  rf <- function(x) { x <- factor(x); relevel(x, ref = if ("CWP&MWA" %in% levels(x)) "CWP&MWA" else "cooler") }
  fitc <- function(f, d) tryCatch({ m <- gam(f, data = d, family = betar(link = "logit"), method = "REML")
    pt <- summary(m)$p.table; pt[grepl("^Thermal", rownames(pt)), c(1, 2, 4), drop = FALSE] }, error = function(e) NULL)
  dplyr::bind_rows(lapply(offsets, function(o) {
    ser <- series_at(o)
    dplyr::bind_rows(lapply(panels, function(p) {
      d  <- ser %>% dplyr::filter(panel == p, !is.na(counts_norm_sv), age > 0.2 + 1e-9)
      dn <- d; dn$Thermal <- rf(as.character(dn$Thermal))
      dl <- d; dl$Thermal <- rf(dl$Thermal_lag)
      specs <- list(A_sage = list(counts_norm_sv ~ s(age) + Thermal, dn),
                    B_sage_k4 = list(counts_norm_sv ~ s(age, k = 4) + Thermal, dn),
                    C_no_sage = list(counts_norm_sv ~ Thermal, dn),
                    D_lag100 = list(counts_norm_sv ~ s(age) + Thermal, dl))
      dplyr::bind_rows(lapply(names(specs), function(s) { r <- fitc(specs[[s]][[1]], specs[[s]][[2]]); if (is.null(r)) return(NULL)
        data.frame(o = o, panel = p, spec = s, level = sub("Thermal", "", rownames(r)), est = r[, 1], se = r[, 2], p = r[, 3]) }))
    }))
  }))
}

## ---- Monte Carlo and bootstrap over all combined candidates (M1-M7) ---------------------
run_mc_full <- function(reps, seed) {
  rec0  <- recs %>% dplyr::select(spot, source, t_ka)
  strat <- paste(rec0$spot, ifelse(rec0$spot == "Florida", rec0$source, ""))
  idx   <- split(seq_len(nrow(rec0)), strat)
  sel_main <- c("Punta Maroma" = best_pm_cb$best_model, "Florida" = best_fl_cb$best_model,
                "Belize" = best_bz_cb$best_model, "Saint Croix" = best_sc_cb$best_model,
                "Regional" = best_reg_cb$best_model)
  one <- function(b, sc) {
    r <- rec0
    if (sc %in% c("Age", "Both")) r$t_ka <- age_draws[, b]
    if (sc %in% c("Sampling", "Both")) r <- r[unlist(lapply(idx, function(i) i[sample.int(length(i), length(i), replace = TRUE)])), ]
    ser <- build_series(r)
    dplyr::bind_rows(lapply(panel_names, function(p) {
      d <- ser %>% dplyr::filter(panel == p, !is.na(counts_norm_sv), !is.na(GMHT_CPS0_30N), !is.na(Hur100_av))
      fits <- Filter(Negate(is.null), lapply(drop_hydro(formulas_cb, d), safe_fit, d = d, combined = TRUE))
      aic <- sapply(fits, `[[`, "aic"); sel <- sel_main[[p]]
      data.frame(Panel = p, Scenario = sc, rep = b, best = names(which.min(aic)),
        dAIC_M1 = if ("M1" %in% names(aic)) unname(aic["M1"] - min(aic)) else NA_real_, sel = sel,
        dAIC_sel = if (!is.null(fits[[sel]])) unname(aic[sel] - min(aic)) else NA_real_,
        p_Thermal = pget(fits, sel, "^Thermal$"), p_Hurricane = pget(fits, sel, "^Hurricane$"),
        p_temp = pget(fits, sel, "GMHT"), p_dt = pget(fits, sel, "dT_raw"), p_hur = pget(fits, sel, "Hur100"))
    }))
  }
  set.seed(seed)
  dplyr::bind_rows(lapply(c("Age", "Sampling", "Both"), function(sc)
    dplyr::bind_rows(lapply(reps, function(b) tryCatch(one(b, sc), error = function(e) NULL)))))
}
