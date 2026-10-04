suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
PAN <- c("Punta Maroma","Florida","Belize","Saint Croix","Regional")
det <- read.csv("det_results.csv")
v4ref <- c("Punta Maroma"="M5_av", Florida="M4_av", Belize="M4_av", "Saint Croix"="M4_fqcy", Regional="M4_av")
cref <- setNames(det$best[det$scen=="CORR_01"], det$panel[det$scen=="CORR_01"])
scs <- intersect(c("A","B","AB","A02","AZ"), sub("res_(.*)\\.rds","\\1", list.files(pattern="^res_.*rds$")))
S <- list(); TP <- list(); GP <- list()
for (s in scs) { r <- readRDS(paste0("res_", s, ".rds")); f <- r$fit
  S[[s]] <- f %>% group_by(panel) %>% summarise(scen = s, B = n(),
     sel_v4 = mean(best == v4ref[panel[1]], na.rm=TRUE), sel_corr = mean(best == cref[panel[1]], na.rm=TRUE),
     sel_M1 = mean(best == "M1", na.rm=TRUE), top = names(sort(table(best), decreasing=TRUE))[1], top_f = max(table(best))/n(),
     dAIC_med = median(dAIC_M1, na.rm=TRUE), dAIC_lo = quantile(dAIC_M1, .025, na.rm=TRUE), dAIC_hi = quantile(dAIC_M1, .975, na.rm=TRUE),
     p_dAIC_gt2 = mean(dAIC_M1 > 2, na.rm=TRUE), dev_med = median(dev_best, na.rm=TRUE), .groups="drop")
  ## term p-values of the v4 reference model
  TP[[s]] <- bind_rows(lapply(seq_len(nrow(f)), function(i) { p <- f$p_v4ref[[i]]; if (is.list(p)) p <- p[[1]]
      if (length(p) == 1 && is.na(p)) return(NULL); data.frame(scen = s, panel = f$panel[i], term = names(p), p = unname(p)) })) %>%
    filter(term != "s(age)") %>% group_by(scen, panel, term) %>% summarise(frac_p05 = mean(p < 0.05), p_med = median(p), n = n(), .groups = "drop")
  if (!is.null(r$gap)) GP[[s]] <- data.frame(scen = s, as.data.frame(r$gap)) %>% select(-rep) %>% summarise(across(everything(), ~mean(.x == 0)))
}
S <- bind_rows(S); TP <- bind_rows(TP); GP <- bind_rows(GP)
write.csv(S, "mc_selection_summary.csv", row.names=FALSE); write.csv(TP, "mc_terms_v4ref.csv", row.names=FALSE); write.csv(GP, "mc_gaps_empty.csv", row.names=FALSE)
options(width = 200); print(as.data.frame(S %>% mutate(across(where(is.numeric), ~round(.x, 2)))))
print(as.data.frame(TP %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
print(GP)
