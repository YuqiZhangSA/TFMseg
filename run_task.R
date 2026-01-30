#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1)

task_id <- as.integer(args[1])

dirs <- c("Method_funcs", "Detection", "Plot")
for (d in dirs) {
  r_files <- list.files(d, pattern = "\\.R$", full.names = TRUE)
  for (f in r_files) source(f)
}

suppressPackageStartupMessages({
  library(tensorMiss)
  library(rTensor)
  library(RTFA)
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

Ts <- c(400, 800, 1600, 3200)
dim_obs_list <- list(c(10,10,10), c(20,20,20), c(10,20,40), c(10,10,100))

simulation_settings <- unlist(
  lapply(dim_obs_list, function(dobs) {
    lapply(Ts, function(Time) {
      list(Time = Time, dim_obs = dobs, dim_latent = c(3, 3, 3))
    })
  }),
  recursive = FALSE
)

methods <- c("TFMseg", "TFMseg.vec", "FMseg", "LR")
nrep <- 100L

n_settings <- length(simulation_settings)
n_methods <- length(methods)
n_tasks <- n_settings * n_methods * nrep

if (task_id < 1L || task_id > n_tasks) stop("task_id out of range")

rep_id <- ((task_id - 1L) %% nrep) + 1L
tmp <- (task_id - 1L) %/% nrep
method_id <- (tmp %% n_methods) + 1L
setting_id <- (tmp %/% n_methods) + 1L

setting <- simulation_settings[[setting_id]]
method <- methods[method_id]

Time <- setting$Time
dim_obs <- setting$dim_obs
dim_latent <- setting$dim_latent

seed_start <- 900L
data_setting <- "s2"
dist <- "Gaussian"
dep <- TRUE
theta_coef <- c(0.25, 0.5, 0.75)
r_hat <- c(3, 3, 3)
trim_coef <- 1/4
threshold_coef <- c(0.11684, 0.70480, 1.12316)
acc_coef <- 1
lrv <- TRUE
V.diag <- TRUE

set.seed(seed_start + rep_id)

res <- simu_comparison(
  method = method,
  nrep = 1,
  seed_start = seed_start + rep_id - 1L,
  data_setting = data_setting,
  dim_obs = dim_obs,
  dim_latent = dim_latent,
  Time = Time,
  dist = dist,
  coeff = 0.7,
  dep = dep,
  m = NULL,
  theta_coef = theta_coef,
  r_hat = r_hat,
  trim_coef = trim_coef,
  threshold_coef = threshold_coef,
  thd.type = if (method %in% c("TFMseg","TFMseg.vec")) {
    if (method == "TFMseg") "fixed" else "oracle"
  } else "fixed",
  V.diag = V.diag,
  lrv = lrv,
  acc_coef = acc_coef
)

dir.create("results_tasks", showWarnings = FALSE, recursive = TRUE)

out_file <- sprintf(
  "results_tasks/task_%05d_method-%s_T-%d_dim-%s_rep-%03d.rds",
  task_id, method, Time, paste(dim_obs, collapse="x"), rep_id
)
saveRDS(list(
  task_id = task_id,
  setting_id = setting_id,
  method = method,
  Time = Time,
  dim_obs = dim_obs,
  rep_id = rep_id,
  res = res
), out_file)

cat("Saved:", out_file, "\n")
