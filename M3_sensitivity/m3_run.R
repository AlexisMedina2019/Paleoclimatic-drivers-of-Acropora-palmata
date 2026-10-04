## M3 sensitivity driver. Usage: Rscript m3_run.R <scenario> <B> <ncores>
args <- commandArgs(TRUE); scen <- args[1]; B <- as.integer(args[2]); nc <- as.integer(args[3])
source("m3_core.R"); suppressMessages(library(parallel))
meta <- read.csv("draws_meta.csv"); DA <- as.matrix(read.csv("draws_ages.csv")); DZ <- as.matrix(read.csv("draws_ages_BZdR.csv"))
GAPS <- list("Punta Maroma" = c(0.71, 1.44), "Belize_1" = c(1.99, 2.72), "Belize_2" = c(3.68, 4.20), "Belize_3" = c(5.54, 6.00))
BWS <- c(0.05, 0.09, 0.15)
REF <- list(v4 = c("Punta Maroma"="M5_av", Florida="M4_av", Belize="M4_av", "Saint Croix"="M4_fqcy", Regional="M4_av"))
cfg <- switch(sub("_.*", "", scen),
  A = list(binw = 0.1, off = 0, ages = "draw", boot = FALSE),
  B = list(binw = 0.1, off = 0, ages = "central", boot = TRUE),
  AB = list(binw = 0.1, off = 0, ages = "draw", boot = TRUE),
  A02 = list(binw = 0.2, off = 0, ages = "draw", boot = FALSE),
  AZ = list(binw = 0.1, off = 0, ages = "drawBZ", boot = FALSE))
set.seed(42)
one <- function(b) {
  t <- switch(cfg$ages, central = meta$t_used_ka, draw = DA[, b], drawBZ = DZ[, b])
  idx <- seq_len(nrow(meta))
  if (cfg$boot) idx <- unlist(lapply(split(idx, paste(meta$site, ifelse(meta$site == "Florida", meta$stratum, ""))),
                                      function(i) i[sample.int(length(i), length(i), replace = TRUE)]))
  recs <- data.frame(site = meta$site[idx], stratum = meta$stratum[idx], t = t[idx])
  P <- build_panels(recs, cfg$binw, cfg$off)
  out <- list(); kd <- list(); ser <- list()
  for (p in names(P)) {
    r <- fit_panel(P[[p]])
    pv <- lapply(REF$v4[p], function(mn) if (!is.null(r$fits[[mn]])) term_p(r$fits[[mn]]) else NA)
    bp <- if (length(r$best)) term_p(r$fits[[r$best]]) else NA
    out[[p]] <- data.frame(rep = b, panel = p, n = nrow(P[[p]]), best = ifelse(length(r$best), r$best, NA),
      dAIC_M1 = unname(r$aic["M1"] - min(r$aic)), dev_best = if (length(r$best)) summary(r$fits[[r$best]])$dev.expl else NA,
      aic = I(list(r$aic)), p_v4ref = I(pv), p_best = I(list(bp)))
    ser[[p]] <- data.frame(rep = b, panel = p, age = P[[p]]$age, y = P[[p]]$counts_norm_sv)
  }
  if (cfg$binw == 0.1) {
    rw <- record_weights(recs, P)
    for (bw in BWS) {
      for (s in SITES) { x <- rw[rw$site == s, ]; kd[[paste(s, bw)]] <- data.frame(rep = b, panel = s, bw = bw, i = 1:601, y = wkde(x$t, x$w_site, bw)) }
      kd[[paste("Reg", bw)]] <- data.frame(rep = b, panel = "Regional", bw = bw, i = 1:601, y = wkde(rw$t, rw$w_reg, bw))
    }
  }
  gap <- sapply(names(GAPS), function(g) { s <- sub("_.*", "", g); sum(recs$site == s & recs$t > GAPS[[g]][1] & recs$t <= GAPS[[g]][2]) })
  list(fit = do.call(rbind, out), ser = do.call(rbind, ser), kde = do.call(rbind, kd), gap = c(rep = b, gap))
}
t0 <- Sys.time()
## rep-level scheduler: nc concurrent forked replicates, each killed after 240 s
res <- vector("list", B); jobs <- list(); nxt <- 1
while (nxt <= B || length(jobs)) {
  while (length(jobs) < nc && nxt <= B) { jobs[[as.character(nxt)]] <- list(j = mcparallel(one(nxt), silent = TRUE), t0 = Sys.time()); nxt <- nxt + 1 }
  Sys.sleep(0.2)
  for (k in names(jobs)) { jb <- jobs[[k]]; x <- mccollect(jb$j, wait = FALSE)
    if (!is.null(x)) { res[[as.integer(k)]] <- if (is.null(x[[1]])) "crashed" else x[[1]]; jobs[[k]] <- NULL; next }
    if (as.numeric(difftime(Sys.time(), jb$t0, units = "secs")) > 240) { tools::pskill(jb$j$pid, tools::SIGKILL); Sys.sleep(0.1); mccollect(jb$j, wait = FALSE); res[[as.integer(k)]] <- "timeout"; jobs[[k]] <- NULL } }
}
bad <- !vapply(res, is.list, TRUE); cat(scen, "lost replicates:", sum(bad), "of", B, "\n")
res <- res[!bad]
saveRDS(list(fit = do.call(rbind, lapply(res, `[[`, "fit")), ser = do.call(rbind, lapply(res, `[[`, "ser")),
             kde = do.call(rbind, lapply(res, `[[`, "kde")), gap = do.call(rbind, lapply(res, `[[`, "gap")), cfg = cfg),
        paste0("res_", scen, ".rds"))
cat(scen, B, "done in", format(Sys.time() - t0), "\n")
