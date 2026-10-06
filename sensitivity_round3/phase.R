## Phase ensemble (pre-specified 2026-10-05): bin origin shifted by o = 0, 10, ..., 90 yr
## (bins in years before 1950 + o). Corals and TC events are re-binned exactly; temperature is
## interpolated linearly at bin centres; regimes take the BP bin containing the bin centre.
## Youngest bin: storms counted from the youngest dated colony (1956 CE, BP -6) back to its lower edge.
suppressMessages({library(dplyr); library(mgcv); library(readxl)})
load("/home/claude/run_mcfull/env.RData")
age_ref <- cli_tab$age   # exact bin ages used by build_series()
s4 <- read_excel("ads5624_table_s4.xlsx", sheet = 1); names(s4)[c(4, 8)] <- c("yva", "EL")
ev_bp <- as.numeric(s4$yva[grepl("^EL", s4$EL)]) - 72             # event ages, years BP
v4 <- read_excel("df_GAM130326_regimes_v4.xlsx", sheet = "df_v1_GAM")
Tbp <- v4$GMHT_CPS0_30N; cen_bp <- seq(50, 5950, 100)
cat_cols <- c("Thermal", "Hydro", "Hurricane", "Hydro_FL", "Hur_FL")
fit_one <- function(f, d, combined) {
  ns <- count_smooth_terms(f); n <- nrow(d)
  if (n < 8 || n <= ns * 3) return(NULL)
  ff <- if (combined) cap_smooth_k(f, safe_k_for(n, ns)) else f
  m <- tryCatch(gam(ff, data = d, family = betar(link = "logit"), method = "REML", select = combined), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  list(aic = AIC(m), p = tryCatch(term_pvalues(m), error = function(e) numeric(0)))
}
safe_fit <- function(f, d, combined) {
  if (any(c("Thermal", "Hurricane", "Hydro") %in% all.vars(f))) {
    j <- parallel::mcparallel(fit_one(f, d, combined), silent = TRUE); r <- parallel::mccollect(j, wait = TRUE)[[1]]
    if (is.null(r) || inherits(r, "try-error")) NULL else r
  } else fit_one(f, d, combined)
}
res <- list()
offs <- as.numeric(strsplit(Sys.getenv("OFFS", "0,10,20,30,40,50,60,70,80,90"), ",")[[1]])
for (o in offs) {
  lo <- 100 * (0:59) - o; hi <- lo + 100                              # bin edges in BP
  fq <- sapply(1:60, function(k) sum(ev_bp >= ifelse(k == 1, -6, lo[k]) & ev_bp < hi[k]))
  av <- sapply(1:60, function(k) sum(ev_bp >= lo[k] - 50 & ev_bp < hi[k] + 50) / 2)   # centred 200-yr window, as in v4
  ctr <- (lo + hi) / 2
  Tn <- approx(cen_bp, Tbp, ctr, rule = 2)$y
  ib <- pmin(60, pmax(1, floor(ctr / 100) + 1))                       # BP bin holding the centre
  cli_tab <<- tibble(age = age_ref, GMHT_CPS0_30N = Tn, dT_raw = c(0, diff(Tn)),
                     Hur100_fqcy = fq, Hur100_av = av, Thermal = v4$Thermal[ib], Hydro_YUC = ifelse(is.na(v4$Hydro[ib]), "no data", v4$Hydro[ib]),
                     Hydro_FL = v4$Hydro_FL[ib], Hur_YUC = v4$Hurricane[ib], Hur_FL = v4$Hur_FL[ib])
  rec <- recs %>% select(spot, source, t_ka) %>% mutate(t_ka = t_ka + o / 1000)
  ser <- build_series(rec)
  for (p in panel_names) {
    d <- ser %>% filter(panel == p, !is.na(counts_norm_sv))
    mg <- lapply(drop_hydro(formulas_mg, d), safe_fit, d = d, combined = FALSE); mg <- mg[!sapply(mg, is.null)]
    cb <- lapply(drop_hydro(formulas_cb, d), safe_fit, d = d, combined = TRUE);  cb <- cb[!sapply(cb, is.null)]
    a_mg <- sapply(mg, `[[`, "aic"); a_cb <- sapply(cb, `[[`, "aic")
    pget <- function(lst, mod, pat) { if (is.null(lst[[mod]])) return(NA_real_); v <- lst[[mod]]$p[grepl(pat, names(lst[[mod]]$p))]; if (length(v)) min(v) else NA_real_ }
    best_cb <- names(which.min(a_cb))
    if (Sys.getenv("DBG") == "1") { cat("\n", p, "fitted:", names(a_cb), "\n"); print(round(sort(a_cb), 2)) }
    res[[length(res) + 1]] <- data.frame(o = o, panel = p, n = nrow(d), best_marg = names(which.min(a_mg)),
      p_hf_marg = pget(mg, "M_hf", "Hur100"), p_hav_marg = pget(mg, "M_hav", "Hur100"),
      dAIC_hav_vs_age = unname(a_mg["M_hav"] - a_mg[which.min(a_mg)]),
      best_comb = best_cb, dAIC_comb_vs_M1 = unname(a_cb["M1"] - min(a_cb)),
      p_hur_best = pget(cb, best_cb, "Hur100"), p_therm_best = pget(cb, best_cb, "^Thermal$"),
      p_hurreg_best = pget(cb, best_cb, "^Hurricane$"),
      p_hf_M4 = pget(cb, "M4_fqcy", "Hur100"), p_hav_M4 = pget(cb, "M4_av", "Hur100"),
      p_therm_M6av = pget(cb, "M6_av", "^Thermal$"))
    cat(o, p, best_cb, "\n")
  }
}
out <- bind_rows(res); write.csv(out, Sys.getenv("OUT", "phase_ensemble.csv"), row.names = FALSE)
