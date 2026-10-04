source("m3_core.R"); meta <- read.csv("draws_meta.csv")
P <- build_panels(data.frame(site = meta$site, stratum = meta$stratum, t = meta$t_used_ka))
for (p in c("Punta Maroma","Regional","Belize")) for (drop in c(FALSE, TRUE)) {
  d <- P[[p]]; if (drop) d <- d[d$age > 0.15, ]
  r <- fit_panel(d, c("M1","M4_av")); tp <- term_p(r$fits$M4_av)
  cat(sprintf("%-12s drop0.1=%-5s n=%d dAIC(M4_av vs M1)=%.1f p(Hur100_av)=%.2g\n", p, drop, nrow(d), r$aic["M1"]-r$aic["M4_av"], tp["s(Hur100_av)"])) }
