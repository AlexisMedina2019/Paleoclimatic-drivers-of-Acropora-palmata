## Round-3 sensitivity analyses (Section 4.7 of regional_analysis_Final.qmd).
## Sourced inside the qmd: uses the objects of the session (recs, age_draws, build_series,
## cli_tab, formulas_mg, formulas_cb, drop_hydro, count_smooth_terms, cap_smooth_k,
## safe_k_for, term_pvalues, panel_names, best_*_cb).

## ---- diagnostic log (2026-10-11): one line per fit, written before fitting, so that if R aborts
## the last line of round3_fit_log.txt names the offset, panel and formula that crashed it.
.r3ctx <- ""
r3_log <- function(...) cat(format(Sys.time(), "%H:%M:%S"), .r3ctx, ..., "\n",
  file = here("sensitivity_round3", "round3_fit_log.txt"), append = TRUE)

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
  ## dev, pt and rank_def added 2026-10-30 (run_frames, run_cps_members); earlier callers use aic and p only
  list(aic = AIC(m), p = tryCatch(term_pvalues(m), error = function(e) numeric(0)),
       dev = summary(m)$dev.expl, pt = summary(m)$p.table[, c(1, 4), drop = FALSE],
       rank_def = m$rank < length(coef(m)))
}
## Windows has no fork. The crash barrier is a persistent background R session (callr): a fit
## with a categorical regime runs there; if mgcv aborts that session or the fit exceeds 60 s,
## the session is replaced and the fit counts as not fitted, as with the forked child on Unix.
## (Added 2026-10-11: on Windows the phase ensemble aborted R at phase 70, Belize, M7_av.)
.r3_rs <- NULL
r3_session <- function() {
  if (!requireNamespace("callr", quietly = TRUE)) stop("install.packages('callr') is required on Windows")
  if (is.null(.r3_rs) || !.r3_rs$is_alive()) {
    .r3_rs <<- callr::r_session$new(wait = TRUE)
    .r3_rs$run(function() suppressPackageStartupMessages(library(mgcv)))
  }
  .r3_rs
}
fit_in_worker <- function(f, d, combined, timeout_s = 60) {
  ns <- count_smooth_terms(f); n <- nrow(d)
  if (n < 8 || n <= ns * 3) return(NULL)
  ff <- if (combined) cap_smooth_k(f, safe_k_for(n, ns)) else f
  environment(ff) <- globalenv()
  rs <- r3_session()
  rs$call(function(ff, d, combined, tp) {
    m <- tryCatch(mgcv::gam(ff, data = d, family = mgcv::betar(link = "logit"), method = "REML", select = combined),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    sm <- summary(m)
    list(aic = stats::AIC(m), p = tryCatch(tp(m), error = function(e) numeric(0)),
         dev = sm$dev.expl, pt = sm$p.table[, c(1, 4), drop = FALSE],
         rank_def = m$rank < length(stats::coef(m)))
  }, args = list(ff, d, combined, term_pvalues))
  if (rs$poll_process(timeout_s * 1000) != "ready") {
    r3_log("TIMEOUT: worker killed, fit counted as not fitted"); rs$kill(); .r3_rs <<- NULL; return(NULL)
  }
  r <- tryCatch(rs$read(), error = function(e) NULL)
  if (is.null(r) || !is.null(r$error)) {
    r3_log("WORKER CRASHED: fit counted as not fitted")
    if (!rs$is_alive()) .r3_rs <<- NULL
    return(NULL)
  }
  r$result
}
safe_fit <- function(f, d, combined) {
  r3_log(if (combined) "cb" else "mg", nrow(d), deparse1(f))
  has_factor <- any(c("Thermal", "Hurricane", "Hydro") %in% all.vars(f))
  if (!has_factor) return(fit_summary(f, d, combined))
  if (.Platform$OS.type != "unix") return(fit_in_worker(f, d, combined))
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
## 2026-10-30 (v5 table, R1.3): the BP table holds the CPS median at the bin centres (50 ... 5950 yr BP)
## plus T_6.05 (median at 6.05 ka). dT is the change from the preceding (older) century into the bin,
## dT = T(centre) - T(centre + 100 yr); approx(rule = 1) returns NA beyond the grid instead of the
## constant edge value that rule = 2 gave (which produced a spurious dT = 0).
ev_bp  <- read.csv(here("sensitivity_round3", "schmitt2025_tableS4a_event_layers.csv"))$varve_years_before_2022 - 72
cli_bp <- cli_tab                  # BP climate table of the main analysis (Section 4.5.2)
age_ref <- cli_tab$age
make_cli <- function(o) {
  lo <- 100 * (0:59) - o; hi <- lo + 100; ctr <- (lo + hi) / 2
  ib <- pmin(60, pmax(1, floor(ctr / 100) + 1))
  T_grid   <- c(seq(50, 5950, 100), 6050)
  T_vals   <- c(cli_bp$GMHT_CPS0_30N, T_6.05)
  Tn       <- approx(T_grid, T_vals, ctr,       rule = 2)$y   # bin centres; youngest shifted centre < 50 yr takes the 50-yr value
  Tn_older <- approx(T_grid, T_vals, ctr + 100, rule = 1)$y
  hy <- as.character(cli_bp$Hydro_YUC)[ib]
  tibble::tibble(age = age_ref, GMHT_CPS0_30N = Tn, dT_raw = Tn - Tn_older,
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
      .r3ctx <<- paste("phase", o, p)
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
      .r3ctx <<- paste("phase-coef", o, p); r3_log("s(age) + Hurricane")
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
      .r3ctx <<- paste("lia", o, p); r3_log("LIA specs")
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
      .r3ctx <<- paste("mc", sc, b, p)
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


## =========================================================================================
## Added 2026-10-30 (round 3, reviewer point R1.3): v5 climate table
## =========================================================================================

## rank-deficient categorical design: excluded before fitting, as in fit_combined_gams_safe()
rank_ok <- function(f, d) {
  fv <- intersect(c("Thermal", "Hydro", "Hurricane"), all.vars(f))
  if (!length(fv)) return(TRUE)
  X <- model.matrix(reformulate(fv), droplevels(d))
  qr(X)$rank == ncol(X)
}
fit_candidates <- function(fset, d, combined = TRUE) {
  fset <- drop_hydro(fset, d)
  fits <- lapply(fset, function(f) if (rank_ok(f, d)) safe_fit(f, d = d, combined = combined) else NULL)
  fits <- Filter(Negate(is.null), fits)
  Filter(function(x) !isTRUE(x$rank_def), fits)
}
fmt_p3 <- function(p) ifelse(p < 1e-4, "<1e-04", sprintf("%.2g", p))

## climate table (same columns as cli_tab) from a climate workbook
cli_from_table <- function(file) {
  x <- as.data.frame(readxl::read_excel(file))
  x$interval <- factor(x$interval, levels = unique(x$interval))
  x %>% dplyr::transmute(age = breaks[match(interval, levels(interval)) + 1], GMHT_CPS0_30N, dT_raw,
                         Hur100_fqcy, Hur100_av, Thermal, Hydro_YUC = ifelse(is.na(Hydro), "no data", Hydro),
                         Hydro_FL, Hur_YUC = Hurricane, Hur_FL)
}

## ---- Section 4.7.1: youngest bin and time reference (three frames, full candidate set) ----
run_frames <- function() {
  frames <- list(
    list(label = "BP, youngest bin to 1956 CE (main)", file = here("df_GAM130326_regimes_v5.xlsx"), shift = 0),
    list(label = "BP, youngest bin to 2022 CE",        file = here("df_GAM130326_regimes_v5c.xlsx"), shift = 0),
    list(label = "Years before 2021 CE",
         file = here("sensitivity_round3", "df_GAM130326_regimes_2021frame_v5.xlsx"), shift = 0.071))
  t2 <- list(); yb <- list()
  for (fr in frames) {
    old <- cli_tab; cli_tab <<- cli_from_table(fr$file)
    ser <- build_series(recs %>% dplyr::select(spot, source, t_ka) %>% dplyr::mutate(t_ka = t_ka + fr$shift))
    cli_tab <<- old
    for (p in panel_names) {
      .r3ctx <<- paste("frame", fr$label, p)
      d    <- ser %>% dplyr::filter(panel == p, !is.na(counts_norm_sv))
      fits <- fit_candidates(formulas_cb, d)
      aic  <- sapply(fits, `[[`, "aic"); sel <- names(which.min(aic))
      w    <- exp(-0.5 * (aic - min(aic))); w <- w / sum(w)
      ps   <- fits[[sel]]$p; ps <- ps[names(ps) != "s(age)"]
      sig  <- ps[ps < 0.05]
      pt   <- fits[[sel]]$pt; pt <- pt[grepl("^(Thermal|Hurricane|Hydro)", rownames(pt)), , drop = FALSE]
      has_fac <- any(c("Thermal", "Hurricane", "Hydro") %in% all.vars(formulas_cb[[sel]]))
      con  <- pt[pt[, 2] < 0.05, , drop = FALSE]
      t2[[length(t2) + 1]] <- data.frame(Frame = fr$label, Panel = p, Selected = sel,
        AIC = round(aic[[sel]], 3), Weight = round(w[[sel]], 2), DevExpl = round(100 * fits[[sel]]$dev, 1),
        Significant = if (length(sig)) paste(sprintf("%s (p=%s)", names(sig), fmt_p3(sig)), collapse = "; ") else "Not significant",
        Contrasts = if (!has_fac) "n/a (no categorical terms)" else if (!nrow(con)) "none" else
          paste(sprintf("%s (%.2f, p=%s)", rownames(con), con[, 1], fmt_p3(con[, 2])), collapse = "; "),
        check.names = FALSE)
      d0 <- d %>% dplyr::filter(age > min(age) + 1e-9)
      f0 <- fit_candidates(formulas_cb[c("M1", sel)], d0)
      yb[[length(yb) + 1]] <- data.frame(Frame = fr$label, Panel = p,
        `dAIC_vs_s(age)_all_bins` = round(aic[["M1"]] - aic[[sel]], 1),
        `dAIC_vs_s(age)_without_youngest` = if (all(c("M1", sel) %in% names(f0)))
          round(f0[["M1"]]$aic - f0[[sel]]$aic, 1) else NA_real_, check.names = FALSE)
    }
  }
  list(table2 = dplyr::bind_rows(t2), youngest = dplyr::bind_rows(yb))
}

## ---- Section 4.7.5: dependence on the CPS ensemble member ---------------------------------
## ens: data frame with binAges (bin centres, yr BP) and members V1 ... V500 (Kaufman et al. 2020,
## CPS 0-30N). For each member the temperature anomaly is interpolated at the bin centres and dT is
## the change from the preceding (older) century, exactly as for the median (Section 1.2). The coral
## series and all other predictors are those of the main analysis (merged_df_g). "median" is
## included as a check: it must reproduce the main analysis.
run_cps_members <- function(ens, n_members = 100, seed = 20261030) {
  age_k <- ens$binAges
  X     <- as.matrix(ens[, grepl("^V[0-9]+$", names(ens))])
  set.seed(seed)
  pick  <- c(1, sample(2:ncol(X), n_members - 1))
  cols  <- c(list(median = apply(X, 1, median)), lapply(setNames(pick, colnames(X)[pick]), function(i) X[, i]))
  ctr   <- seq(50, 5950, 100)
  sel_main <- c("Punta Maroma" = best_pm_cb$best_model, "Florida" = best_fl_cb$best_model,
                "Belize" = best_bz_cb$best_model, "Saint Croix" = best_sc_cb$best_model,
                "Regional" = best_reg_cb$best_model)
  base  <- merged_df_g %>% dplyr::filter(!is.na(counts_norm_sv)) %>% dplyr::mutate(age_key = round(age, 1))
  sage  <- sapply(panel_names, function(p) {
    .r3ctx <<- paste("cps s(age)", p)
    safe_fit(counts_norm_sv ~ s(age), base %>% dplyr::filter(panel == p), combined = FALSE)$aic })
  dplyr::bind_rows(lapply(names(cols), function(mb) {
    v   <- cols[[mb]]
    Tn  <- approx(age_k, v, ctr)$y; To <- approx(age_k, v, ctr + 100)$y
    tab <- data.frame(age_key = round(ctr / 1000 + 0.05, 1), GMHT_CPS0_30N = Tn, dT_raw = Tn - To)
    dplyr::bind_rows(lapply(panel_names, function(p) {
      .r3ctx <<- paste("cps", mb, p)
      d  <- base %>% dplyr::filter(panel == p) %>% dplyr::select(-GMHT_CPS0_30N, -dT_raw) %>%
        dplyr::left_join(tab, by = "age_key")
      mt <- safe_fit(formulas_mg$M_temp, d, combined = FALSE)
      md <- safe_fit(formulas_mg$M_dt,   d, combined = FALSE)
      cb <- fit_candidates(formulas_cb, d)
      aic <- sapply(cb, `[[`, "aic"); sel <- sel_main[[p]]
      data.frame(member = mb, Panel = p,
        dAIC_temp_vs_sage = if (is.null(mt)) NA_real_ else mt$aic - sage[[p]],
        p_temp_marg = if (is.null(mt)) NA_real_ else unname(mt$p["s(GMHT_CPS0_30N)"]),
        dAIC_dt_vs_sage = if (is.null(md)) NA_real_ else md$aic - sage[[p]],
        p_dt_marg = if (is.null(md)) NA_real_ else unname(md$p["s(dT_raw)"]),
        best = names(which.min(aic)), sel_main = sel,
        dAIC_sel_vs_M1 = if (all(c("M1", sel) %in% names(aic))) unname(aic["M1"] - aic[sel]) else NA_real_,
        p_temp_sel = pget(cb, sel, "GMHT"), p_dt_sel = pget(cb, sel, "dT_raw"),
        p_Thermal_sel = pget(cb, sel, "^Thermal$"), p_Hurricane_sel = pget(cb, sel, "^Hurricane$"),
        ## model selected with this member
        best_has_regime = any(c("Thermal", "Hurricane") %in% all.vars(formulas_cb[[names(which.min(aic))]])),
        p_temp_best = pget(cb, names(which.min(aic)), "GMHT"), p_dt_best = pget(cb, names(which.min(aic)), "dT_raw"),
        p_Thermal_best = pget(cb, names(which.min(aic)), "^Thermal$"),
        p_Hurricane_best = pget(cb, names(which.min(aic)), "^Hurricane$"),
        ## AIC distance from the best model to the best candidate containing a regime (0 = a regime model is best)
        dAIC_best_regime = { rg <- names(aic)[vapply(names(aic), function(m)
            any(c("Thermal", "Hurricane") %in% all.vars(formulas_cb[[m]])), logical(1))]
          if (length(rg)) min(aic[rg]) - min(aic) else NA_real_ })
    }))
  }))
}
