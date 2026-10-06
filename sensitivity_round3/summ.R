suppressMessages(library(dplyr))
fs <- Sys.glob(c("/home/claude/run_mcfull/M3_sensitivity/gam_uncertainty_mc_full_v4b_part*.rds","/home/claude/run_mcfull2/M3_sensitivity/gam_uncertainty_mc_full_v4b_part*.rds"))
mc <- bind_rows(lapply(fs, readRDS)); cat("files:", length(fs), " rows:", nrow(mc), "\n")
pn <- c("Punta Maroma","Florida","Belize","Saint Croix","Regional")
tb <- mc %>% group_by(Panel = factor(Panel, pn), Scenario = factor(Scenario, c("Age","Sampling","Both"))) %>%
  summarise(n = n(), sel = first(sel),
    top = { t <- sort(table(best), decreasing = TRUE); sprintf("%s (%.0f%%)", names(t)[1], 100*t[1]/n()) },
    sel_best = round(100*mean(best == sel)), sel_2AIC = round(100*mean(dAIC_sel <= 2, na.rm = TRUE)),
    regime_best = round(100*mean(grepl("M5|M6|M7", best))), gt_age = round(100*mean(dAIC_M1 > 2, na.rm=TRUE)),
    Thermal = round(100*mean(p_Thermal < .05, na.rm = TRUE)), HurReg = round(100*mean(p_Hurricane < .05, na.rm = TRUE)),
    GMHT = round(100*mean(p_temp < .05, na.rm = TRUE)), dT = round(100*mean(p_dt < .05, na.rm = TRUE)),
    Hur100 = round(100*mean(p_hur < .05, na.rm = TRUE)), .groups = "drop")
options(width = 200); print(as.data.frame(tb))
write.csv(tb, "/home/claude/run_mcfull/mc_full_summary.csv", row.names = FALSE)
