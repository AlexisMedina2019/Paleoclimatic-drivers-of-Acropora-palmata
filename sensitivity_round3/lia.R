## Pre-specified LIA tests (2026-10-05, 22:00): run once, reported as they come.
suppressMessages({library(dplyr); library(mgcv); library(readxl)})
load("/home/claude/run_mcfull/env.RData"); age_ref <- cli_tab$age
s4 <- read_excel("/home/claude/run_phase/ads5624_table_s4.xlsx", sheet = 1); names(s4)[c(4, 8)] <- c("yva", "EL")
ev_bp <- as.numeric(s4$yva[grepl("^EL", s4$EL)]) - 72
v4 <- read_excel("/home/claude/run_phase/df_GAM130326_regimes_v4.xlsx", sheet = "df_v1_GAM")
cen_bp <- seq(50, 5950, 100)
mk_cli <- function(o) {
  lo <- 100 * (0:59) - o; hi <- lo + 100; ctr <- (lo + hi) / 2
  ib <- pmin(60, pmax(1, floor(ctr / 100) + 1)); Tn <- approx(cen_bp, v4$GMHT_CPS0_30N, ctr, rule = 2)$y
  tibble(age = age_ref, GMHT_CPS0_30N = Tn, dT_raw = c(0, diff(Tn)),
         Hur100_fqcy = sapply(1:60, function(k) sum(ev_bp >= ifelse(k == 1, -6, lo[k]) & ev_bp < hi[k])),
         Hur100_av = sapply(1:60, function(k) sum(ev_bp >= lo[k] - 50 & ev_bp < hi[k] + 50) / 2),
         Thermal = v4$Thermal[ib], Hydro_YUC = ifelse(is.na(v4$Hydro[ib]), "no data", v4$Hydro[ib]),
         Hydro_FL = v4$Hydro_FL[ib], Hur_YUC = v4$Hurricane[ib], Hur_FL = v4$Hur_FL[ib],
         Thermal_lag = c(v4$Thermal[ib][-1], tail(v4$Thermal[ib], 1)))   # regime of the next-older bin (100-yr lag)
}
fitc <- function(f, d) tryCatch({ m <- gam(f, data = d, family = betar(link = "logit"), method = "REML")
  pt <- summary(m)$p.table; pt[grepl("^Thermal", rownames(pt)), c(1, 2, 4), drop = FALSE] }, error = function(e) NULL)
res <- list(); raw <- list()
for (o in seq(0, 90, 10)) {
  cli_tab <- mk_cli(o)
  ser <- build_series(recs %>% select(spot, source, t_ka) %>% mutate(t_ka = t_ka + o / 1000))
  for (p in c("Punta Maroma", "Belize", "Regional")) {
    d <- ser %>% filter(panel == p, !is.na(counts_norm_sv)) %>% mutate(Thermal = as.character(Thermal))
    if (o == 0) raw[[p]] <- d %>% group_by(Thermal) %>% summarise(bins = n(), mean_y = mean(counts_norm_sv), median_y = median(counts_norm_sv), .groups = "drop") %>% mutate(panel = p)
    dn <- d %>% filter(age > 0.2 + 1e-9)                                # without the two modern centuries
    rf <- function(x) { x <- factor(x); relevel(x, ref = if ("CWP&MWA" %in% levels(x)) "CWP&MWA" else "cooler") }
    dn$Thermal <- rf(dn$Thermal)
    dl <- d %>% filter(age > 0.2 + 1e-9); dl$Thermal <- rf(dl$Thermal_lag)
    specs <- list(A_sage = list(counts_norm_sv ~ s(age) + Thermal, dn),
                  B_sage_k4 = list(counts_norm_sv ~ s(age, k = 4) + Thermal, dn),
                  C_no_sage = list(counts_norm_sv ~ Thermal, dn),
                  D_lag100 = list(counts_norm_sv ~ s(age) + Thermal, dl))
    for (s in names(specs)) { r <- fitc(specs[[s]][[1]], specs[[s]][[2]]); if (is.null(r)) next
      res[[length(res) + 1]] <- data.frame(o = o, panel = p, spec = s, level = sub("Thermal", "", rownames(r)), est = r[, 1], se = r[, 2], p = r[, 3]) }
  }
  cat(o, "\n")
}
write.csv(bind_rows(raw), "raw_by_regime.csv", row.names = FALSE); write.csv(bind_rows(res), "lia_tests.csv", row.names = FALSE)
