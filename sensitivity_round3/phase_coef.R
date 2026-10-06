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
  for (p in c("Punta Maroma", "Belize", "Regional")) {
    d <- ser %>% filter(panel == p, !is.na(counts_norm_sv))
    d$Hurricane <- relevel(droplevels(d$Hurricane), ref = "moderate")
    j <- parallel::mcparallel({ m <- gam(counts_norm_sv ~ s(age) + Hurricane, data = d, family = betar(link = "logit"), method = "REML")
      cf <- summary(m)$p.table; cf[grepl("^Hurricane", rownames(cf)), c(1, 4), drop = FALSE] }, silent = TRUE)
    r <- parallel::mccollect(j, wait = TRUE)[[1]]
    if (!is.null(r) && !inherits(r, "try-error"))
      res[[length(res) + 1]] <- data.frame(o = o, panel = p, level = sub("Hurricane", "", rownames(r)), est = r[, 1], p = r[, 2])
    tb <- table(d$Hurricane); res[[length(res)]]$n_bins <- as.integer(tb[res[[length(res)]]$level])
  }
}
out <- bind_rows(res); write.csv(out, "phase_regime_coefs.csv", row.names = FALSE)
