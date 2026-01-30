#!/usr/bin/env Rscript

## ---- set project root (assume you run this from ~/TFMseg)
proj_root <- getwd()

dirs <- file.path(proj_root, c("Method_funcs", "Detection", "Plot"))
for (d in dirs) {
  r_files <- list.files(d, pattern = "\\.R$", full.names = TRUE)
  for (f in r_files) source(f)
}

suppressPackageStartupMessages({
  library(tensorMiss)
  library(rTensor)
  library(RTFA)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(kableExtra)
  library(lubridate)
})

Ts <- c(400, 800, 1600, 3200, 6400)

dim_obs_list <- list(
  c(10, 10, 10),
  c(20, 20, 20),
  c(10, 20, 40),
  c(10, 10, 100)
)

simulation_settings <- unlist(
  lapply(dim_obs_list, function(dobs) {
    lapply(Ts, function(Time) {
      list(Time = Time, dim_obs = dobs, dim_latent = c(3, 3, 3))
    })
  }),
  recursive = FALSE
)

res_true_r <- simu_ret(
  methods = c("TFMseg", "TFMseg.vec", "FMseg", "LR"),
  thd.type_values = c("fixed"),
  V_shap_values = c("diag"),
  lrv = TRUE,
  simu_set_detect = simulation_settings,
  dist = "Gaussian",
  data_setting = "s2",
  dep = TRUE,
  theta_coef = c(0.25, 0.5, 0.75),
  r_hat = c(3, 3, 3),
  trim_coef = 1/4,
  nrep = 100,
  threshold_coef = c(0.11684, 0.70480, 1.12316),
  acc_coef = 1
)

## Save outputs
dir.create(file.path(proj_root, "results"), showWarnings = FALSE, recursive = TRUE)

saveRDS(res_true_r, file = file.path(proj_root, "results", "res_true_r.rds"))
write.csv(res_true_r$table, file = file.path(proj_root, "results", "res_true_r_table.csv"), row.names = FALSE)

print(res_true_r$table)
