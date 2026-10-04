source("m3_core.R"); suppressMessages({library(ggplot2); library(tidyr); library(patchwork)})
PAN <- c("Punta Maroma","Florida","Belize","Saint Croix","Regional")
meta <- read.csv("draws_meta.csv"); CORR <- data.frame(site = meta$site, stratum = meta$stratum, t = meta$t_used_ka)
ORIG <- read.csv("records_v4_original.csv")
dser <- readRDS("det_series.rds")
A <- readRDS("res_A.rds"); B <- readRDS("res_B.rds")
env <- function(r, lab) { s <- r$ser; reps <- unique(s$rep)
  s %>% tidyr::complete(nesting(rep), panel, age = round(seq(0.1, 6, 0.1), 1), fill = list(y = 0)) %>%
    group_by(panel, age) %>% summarise(lo = quantile(y, .025), q25 = quantile(y, .25), med = median(y), q75 = quantile(y, .75), hi = quantile(y, .975), .groups = "drop") %>% mutate(src = lab) }
E <- bind_rows(env(A, "Age uncertainty (Monte Carlo)"), env(B, "Sampling uncertainty (bootstrap)"))
pts <- dser %>% filter(scen %in% c("CORR_01", "ORIG_01")) %>% mutate(scen = recode(scen, CORR_01 = "Corrected ages", ORIG_01 = "Submitted ages"))
E$panel <- factor(E$panel, PAN); pts$panel <- factor(pts$panel, PAN)
f1 <- ggplot(E, aes(age)) + geom_ribbon(aes(ymin = lo, ymax = hi, fill = src), alpha = .35) +
  geom_point(data = pts, aes(age, y, shape = scen, colour = scen), size = 1.1) +
  scale_shape_manual(values = c(16, 4)) + scale_colour_manual(values = c("black", "#D7301F")) +
  scale_fill_manual(values = c("#2171B5", "#FDAE6B")) + scale_x_reverse(breaks = 0:6) +
  facet_grid(panel ~ src) + labs(x = "Age (ka BP)", y = "counts_norm_sv", fill = NULL, shape = NULL, colour = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave("FigS_M3_series_envelopes.png", f1, width = 9, height = 9, dpi = 200)
## KDE
grid <- seq(0, 6, by = 0.01)
P0 <- build_panels(CORR); rw <- record_weights(CORR, P0)
cen <- bind_rows(lapply(c(0.05, 0.09, 0.15), function(bw) bind_rows(c(lapply(SITES, function(s) { x <- rw[rw$site == s, ]; data.frame(panel = s, bw = bw, t = grid, y = wkde(x$t, x$w_site, bw)) }),
   list(data.frame(panel = "Regional", bw = bw, t = grid, y = wkde(rw$t, rw$w_reg, bw)))))))
kenv <- function(r, lab) r$kde %>% group_by(panel, bw, i) %>% summarise(lo = quantile(y, .025), hi = quantile(y, .975), .groups = "drop") %>% mutate(t = grid[i], src = lab)
K <- bind_rows(kenv(A, "Age uncertainty"), kenv(B, "Sampling uncertainty"))
K$panel <- factor(K$panel, PAN); cen$panel <- factor(cen$panel, PAN)
gaps <- data.frame(panel = factor(c("Punta Maroma","Belize","Belize","Belize","Regional","Regional"), PAN), x1 = c(.638,1.990,3.684,5.541,2.20,5.10), x2 = c(1.382,2.879,4.202,6.00,3.00,5.90))
f2 <- ggplot() + geom_rect(data = gaps, aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf), fill = "grey85") +
  geom_ribbon(data = K, aes(t, ymin = lo, ymax = hi, fill = src), alpha = .35) + geom_line(data = cen, aes(t, y), linewidth = .4) +
  scale_fill_manual(values = c("#2171B5", "#FDAE6B")) + scale_x_reverse(breaks = 0:6) +
  facet_grid(panel ~ paste("bw =", bw, "ka"), scales = "free_y") + labs(x = "Age (ka BP)", y = "Weighted density", fill = "95 % envelope") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave("FigS_M3_kde_bandwidths.png", f2, width = 10, height = 9, dpi = 200)
cat("figs ok\n")
