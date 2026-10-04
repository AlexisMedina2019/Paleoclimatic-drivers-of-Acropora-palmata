## M3 core: build the modelled panels from a record table exactly as the qmd does, for any bin width,
## and fit the published candidate set. Used by the age-uncertainty Monte Carlo and the bin-width tests.
suppressMessages({library(mgcv); library(dplyr); library(readxl)})
SITES <- c("Punta Maroma", "Florida", "Belize", "Saint Croix")
cli0 <- as.data.frame(read_excel("df_GAM130326_regimes_v3.xlsx"))
cli0$age01 <- round(seq(0.1, 6, 0.1), 1)
HUR_LEV <- c("low", "moderate", "calm", "enhanced", "high", "quiet", "active")

clim_for_bins <- function(binw, offset = 0, shift = 0) {
  cli0 <- cli0; cli0$age01 <- round(cli0$age01 + shift, 1)
  ## climate per bin: numeric = mean of the 0.1-ka rows inside the bin; regimes = row at the bin midpoint
  ub <- seq(binw + offset, 6 + 1e-9, by = binw); lb <- ub - binw
  out <- lapply(seq_along(ub), function(i) {
    rr <- cli0[cli0$age01 > lb[i] + 1e-9 & cli0$age01 <= ub[i] + 1e-9, ]
    if (!nrow(rr)) return(NULL)
    mid <- (lb[i] + ub[i]) / 2; rm <- cli0[which.min(abs((cli0$age01 - 0.05) - mid)), ]
    data.frame(age = round(ub[i], 3), GMHT_CPS0_30N = mean(rr$GMHT_CPS0_30N), dT_raw = mean(rr$dT_raw),
               Hur100_fqcy = mean(rr$Hur100_fqcy), Hur100_av = mean(rr$Hur100_av),
               Thermal = rm$Thermal, Hydro_YUC = ifelse(is.na(rm$Hydro), "no data", rm$Hydro),
               Hydro_FL = rm$Hydro_FL, Hurricane = rm$Hurricane, Hur_FL = rm$Hur_FL)
  })
  bind_rows(out)
}

sv <- function(x) { n <- length(x); y <- round(x / max(x), 4); ifelse(y %in% c(0, 1), (y * (n - 1) + 0.5) / n, y) }

build_panels <- function(recs, binw = 0.1, offset = 0, cli = NULL) {
  ## recs: data.frame(site, stratum, t) with t in ka BP
  if (is.null(cli)) cli <- clim_for_bins(binw, offset)
  ub <- round(seq(binw + offset, 6 + 1e-9, by = binw), 3)
  recs <- recs %>% filter(t <= 6, t > -0.1) %>%
    mutate(age = ub[pmin(length(ub), pmax(1, ceiling((pmax(t, 1e-6) - offset) / binw - 1e-9)))])
  site_series <- function(s) {
    d <- recs %>% filter(site == s)
    if (s == "Florida") {
      d %>% count(stratum, age) %>% group_by(stratum) %>% mutate(std = sv(n)) %>% ungroup() %>%
        group_by(age) %>% summarise(counts = sum(std), .groups = "drop")
    } else d %>% count(age, name = "counts")
  }
  recs$w <- NA_real_
  panels <- lapply(SITES, function(s) { x <- site_series(s); x$counts_norm_sv <- sv(x$counts); x$panel <- s; x })
  names(panels) <- SITES
  reg <- bind_rows(panels) %>% group_by(age) %>% summarise(counts = sum(counts_norm_sv), .groups = "drop")
  reg$counts_norm_sv <- sv(reg$counts); reg$panel <- "Regional"; panels$Regional <- reg
  lapply(panels, function(p) {
    p <- left_join(p, cli, by = "age")
    fl <- p$panel[1] == "Florida"
    p$Hydro <- factor(if (fl) p$Hydro_FL else p$Hydro_YUC, levels = c("wet", "dry", "no data"))
    p$Hurricane <- factor(if (fl) p$Hur_FL else p$Hurricane, levels = HUR_LEV)
    p$Thermal <- relevel(factor(p$Thermal), ref = "warmer")
    p
  })
}

formulas_cb <- list(
  M1 = counts_norm_sv ~ s(age),
  M2 = counts_norm_sv ~ s(age) + s(GMHT_CPS0_30N),
  M3 = counts_norm_sv ~ s(age) + s(GMHT_CPS0_30N) + s(dT_raw),
  M4_fqcy = counts_norm_sv ~ s(age) + s(GMHT_CPS0_30N) + s(dT_raw) + s(Hur100_fqcy),
  M5_fqcy = counts_norm_sv ~ s(age) + s(GMHT_CPS0_30N) + s(dT_raw) + s(Hur100_fqcy) + Hydro,
  M6_fqcy = counts_norm_sv ~ s(age) + s(GMHT_CPS0_30N) + s(dT_raw) + s(Hur100_fqcy) + Hydro + Thermal,
  M7_fqcy = counts_norm_sv ~ s(age) + s(GMHT_CPS0_30N) + s(dT_raw) + s(Hur100_fqcy) + Hydro + Thermal + Hurricane,
  M4_av = counts_norm_sv ~ s(age) + s(GMHT_CPS0_30N) + s(dT_raw) + s(Hur100_av),
  M5_av = counts_norm_sv ~ s(age) + s(GMHT_CPS0_30N) + s(dT_raw) + s(Hur100_av) + Hydro,
  M6_av = counts_norm_sv ~ s(age) + s(GMHT_CPS0_30N) + s(dT_raw) + s(Hur100_av) + Hydro + Thermal,
  M7_av = counts_norm_sv ~ s(age) + s(GMHT_CPS0_30N) + s(dT_raw) + s(Hur100_av) + Hydro + Thermal + Hurricane)

safe_k_for <- function(n_obs, n_smooth) max(3, min(10, floor((n_obs - 1) / max(1, n_smooth))))
cap_k <- function(f, k) as.formula(gsub("s\\(([A-Za-z0-9_]+)\\)", paste0("s(\\1, k = ", k, ")"), deparse(f, width.cutoff = 500)))

fit_one <- function(f, d) {
  ns <- length(gregexpr("s\\(", deparse(f, width.cutoff = 500))[[1]]); n <- nrow(d)
  fv <- intersect(c("Thermal", "Hydro", "Hurricane"), all.vars(f))
  npar <- 1
  if (length(fv)) {
    X <- model.matrix(reformulate(fv), droplevels(d)); if (qr(X)$rank < ncol(X)) return(NULL); npar <- ncol(X)
  }
  if (n < 8 || n <= ns * 3) return(NULL)
  k <- safe_k_for(n, ns)
  tryCatch(if (n > 100) bam(cap_k(f, k), data = d, family = betar(link = "logit"), method = "fREML", select = TRUE, discrete = TRUE)
           else gam(cap_k(f, k), data = d, family = betar(link = "logit"), method = "REML", select = TRUE),
           error = function(e) NULL)
}

term_p <- function(m) {
  s <- summary(m); p <- c(setNames(s$s.table[, "p-value"], rownames(s$s.table)))
  if (!is.null(s$pTerms.table)) p <- c(p, setNames(s$pTerms.table[, "p-value"], rownames(s$pTerms.table)))
  p
}

## run an expression in a forked child; NULL if it crashes, errors or exceeds `secs` (child killed)
run_timeout <- function(expr, secs) {
  j <- parallel::mcparallel(expr, silent = TRUE); t0 <- Sys.time(); r <- NULL
  repeat {
    x <- parallel::mccollect(j, wait = FALSE, timeout = 0.2)
    if (!is.null(x)) { r <- x[[1]]; break }
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > secs) { tools::pskill(j$pid, tools::SIGKILL); Sys.sleep(0.1); parallel::mccollect(j, wait = FALSE); break }
  }
  if (is.null(r) || inherits(r, "try-error")) NULL else r
}
MV <- c("counts_norm_sv","age","GMHT_CPS0_30N","dT_raw","Hur100_fqcy","Hur100_av","Hydro","Thermal","Hurricane")
fit_panel <- function(d, which = names(formulas_cb)) {
  d <- as.data.frame(d); d <- d[complete.cases(d[, MV]), ]
  risky <- nrow(d) < 30 || length(unique(d$counts_norm_sv)) <= 4
  fits <- lapply(formulas_cb[which], function(f) {
    if (!risky) return(fit_one(f, d))
    run_timeout(fit_one(f, d), 60) }); fits <- fits[!vapply(fits, is.null, TRUE)]
  fits <- fits[vapply(fits, function(m) m$rank >= length(coef(m)), TRUE)]
  aic <- vapply(fits, AIC, 1)
  list(fits = fits, aic = aic, best = names(which.min(aic)))
}

## record-level KDE weights, as in qmd §1.4: weights of the records in a bin add up to counts_norm_sv of that bin
record_weights <- function(recs, P, binw = 0.1, offset = 0) {
  ub <- P$Regional$age; ub <- sort(unique(c(ub, unlist(lapply(P[SITES], `[[`, "age")))))
  allub <- seq(binw + offset, 6 + 1e-9, by = binw)
  recs <- recs %>% filter(t <= 6, t > -0.1) %>%
    mutate(age = round(allub[pmin(length(allub), pmax(1, ceiling((pmax(t, 1e-6) - offset) / binw - 1e-9)))], 3))
  fl <- recs %>% filter(site == "Florida") %>% count(stratum, age) %>% group_by(stratum) %>% mutate(fstd = sv(n) / n) %>% ungroup()
  recs <- recs %>% left_join(fl %>% select(stratum, age, fstd), by = c("stratum", "age")) %>%
    left_join(bind_rows(lapply(SITES, function(s) data.frame(site = s, age = P[[s]]$age, fsite = P[[s]]$counts_norm_sv / P[[s]]$counts))), by = c("site", "age")) %>%
    left_join(data.frame(age = P$Regional$age, freg = P$Regional$counts_norm_sv / P$Regional$counts), by = "age") %>%
    mutate(w_site = ifelse(site == "Florida", fstd, 1) * fsite, w_reg = w_site * freg)
  recs
}

wkde <- function(t, w, bw, grid = seq(0, 6, by = 0.01)) {
  d <- density(t, weights = w / sum(w), bw = bw, from = 0, to = 6, n = length(grid)); d$y
}
