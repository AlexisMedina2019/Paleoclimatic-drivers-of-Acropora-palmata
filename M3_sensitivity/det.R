source("m3_core.R")
e <- new.env(); load("../../v4/final.RData", envir = e); r0 <- as.data.frame(e$centennial_counts_sites)
ORIG <- data.frame(site = r0$spot, stratum = r0$source, t = r0$t_ka); write.csv(ORIG, "records_v4_original.csv", row.names = FALSE)
meta <- read.csv("draws_meta.csv")
CORR <- data.frame(site = meta$site, stratum = meta$stratum, t = meta$t_used_ka)
CZ <- CORR; bz <- meta$site == "Belize" & meta$method == "14C"; CZ$t[bz] <- meta$t_sens_dR_regional_ka[bz]
sc <- list(ORIG_01 = list(ORIG, 0.1, 0, 0), CORR_01 = list(CORR, 0.1, 0, 0), CORR_BZdR = list(CZ, 0.1, 0, 0),
           ORIG_02a = list(ORIG, 0.2, 0, 0), CORR_02a = list(CORR, 0.2, 0, 0), CORR_02b = list(CORR, 0.2, 0.1, 0),
           CORR_climP = list(CORR, 0.1, 0, 0.1), CORR_climM = list(CORR, 0.1, 0, -0.1))
out <- list(); ser <- list()
for (s in names(sc)) { a <- sc[[s]]
  P <- build_panels(a[[1]], a[[2]], a[[3]], cli = clim_for_bins(a[[2]], a[[3]], a[[4]]))
  for (p in names(P)) { r <- fit_panel(P[[p]]); b <- r$best; w <- exp(-0.5*(r$aic - min(r$aic))); w <- w/sum(w)
    tp <- term_p(r$fits[[b]]); tp <- tp[names(tp) != "s(age)"]
    out[[paste(s,p)]] <- data.frame(scen = s, panel = p, n = sum(complete.cases(P[[p]][, MV])), best = b, w_best = round(max(w),2),
      dAIC_M1 = round(r$aic["M1"] - min(r$aic), 1), dev = round(summary(r$fits[[b]])$dev.expl*100, 1),
      climate_terms_p = paste(sprintf("%s=%.3g", names(tp), tp), collapse = "; "),
      aic_all = paste(sprintf("%s:%.1f", names(r$aic), r$aic), collapse=" "))
    ser[[paste(s,p)]] <- data.frame(scen = s, panel = p, age = P[[p]]$age, y = P[[p]]$counts_norm_sv)
  } }
D <- do.call(rbind, out); write.csv(D, "det_results.csv", row.names = FALSE); saveRDS(do.call(rbind, ser), "det_series.rds")
print(D[, c("scen","panel","n","best","w_best","dAIC_M1","dev","climate_terms_p")], row.names = FALSE, right = FALSE)
