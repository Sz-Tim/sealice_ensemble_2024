# Compile output from long runs
# Sea lice ms 2024
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk

# These simulations are intended for visualizations of 2023 vertical distributions.



# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(furrr)
library(carrier)
library(sf)
library(rstan)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
theme_set(theme_bw() + theme(panel.grid=element_blank()))

mesh_i <- st_read("data/WeStCOMS2_mesh.gpkg") |>
  mutate(vol_top50m=area*pmin(depth, 50)) |> 
  st_drop_geometry() |> 
  select(i, area, vol_top50m)
out_dir <- "out/sim_2023/"
sim_i <- read_csv(glue("{out_dir}/sim_i.csv")) |>
  mutate(sim=paste0("sim_", i)) |>
  filter(!fixDepth)
site_i <- read_csv("data/farm_sites_2023.csv")

ncores <- 25




# ensemble weights --------------------------------------------------------

mod <- "ranef"
out_ensMixRE <- readRDS(glue("out/ensembles/ensMix_3D_{mod}_stanfit.rds"))
dat_ensMixRE <- readRDS(glue("out/ensembles/ensMix_3D_{mod}_standata.rds"))

iter_sub <- sample.int(length(rstan::extract(out_ensMixRE, pars="sigma")[[1]]), size=3000)
b_p_post <- rstan::extract(out_ensMixRE, pars="b_p")[[1]]

rm(out_ensMixRE); gc()





# calculate ensemble distributions ----------------------------------------

dates <- seq(ymd("2023-01-01"), ymd("2023-12-31"), by=1)
for(i in seq_along(dates)) {
  
  if(file.exists(glue("{out_dir}/processed/verticalDistributions_{dates[i]}.rds"))) next
  cat("\nStarting", i, "of", length(dates), "\n")
  z_df <- load_vertDistr_simSets(out_dir, mesh_i |> select(i),
                                 sim_i,
                                 liceScale=1,
                                 ncores=ncores,
                                 stage="Mature",
                                 date_rng=rep(dates[i], 2)) |>
    pivot_wider(names_from="sim", values_from="value") 
  
  if(nrow(z_df) == 0) next
  
  cat("\nEnsembling", i, "of", length(dates), "\n")
  plan(multicore, workers=ncores)
  z_mx <- z_df |> select(all_of(dat_ensMixRE$sim_names)) |> as.matrix()
  z_mx[is.na(z_mx)] <- 0
  z_df <- z_df |>
    mutate(row=row_number()) |>
    mutate(ens=future_map(row, ~apply(b_p_post[iter_sub,], 1, function(x) sum(x * z_mx[.x,]))),
           ens_mn=map_dbl(ens, mean),
           ens_sd=map_dbl(ens, sd)) |>
    select(i, z, hour, starts_with("ens_"))
  plan(sequential)
  saveRDS(z_df, glue("{out_dir}/processed/verticalDistributions_{dates[i]}.rds"))
  gc()
}





# summarise by region -----------------------------------------------------

f <- dirf(glue("{out_dir}lprocessed"), "verticalDistributions_")
z_ls <- z_linnhe <- z_skye <- vector("list", 12)
for(m in 1:12) {
  z_df <- readRDS(glue("{out_dir}/processed_2023/verticalDistributions_{month.abb[m]}.rds"))
  gc()
  
  z_linnhe[[m]] <- z_df |>
    filter(i %in% linnhe_i$i) |>
    mutate(time=ymd_hms("2023-01-01 00:00:00") + dhours(as.numeric(hour) - 1),
           day=date(time)) |>
    group_by(sim, day, z) |>
    summarise(N=sum(value)/24) |>
    group_by(sim, day) |>
    mutate(prop=N/sum(N))
  gc()
  
  z_skye[[m]] <- z_df |>
    filter(i %in% skye_mesh$i) |>
    mutate(time=ymd_hms("2023-01-01 00:00:00") + dhours(as.numeric(hour) - 1),
           day=date(time)) |>
    group_by(sim, day, z) |>
    summarise(N=sum(value)/24) |>
    group_by(sim, day) |>
    mutate(prop=N/sum(N))
  gc()
  
  z_ls[[m]] <- z_df |>
    mutate(time=ymd_hms("2023-01-01 00:00:00") + dhours(as.numeric(hour) - 1),
           day=date(time)) |>
    group_by(sim, day, z) |>
    summarise(N=sum(value)/24) |>
    group_by(sim, day) |>
    mutate(prop=N/sum(N))
  rm(z_df)
  gc()
}
z_linnhe |>
  reduce(bind_rows) |>
  saveRDS(glue("{out_dir}/processed/summary_daily_z_Linnhe.rds"))
z_skye |>
  reduce(bind_rows) |>
  saveRDS(glue("{out_dir}/processed/summary_daily_z_Skye.rds"))
z_ls |>
  reduce(bind_rows) |>
  saveRDS(glue("{out_dir}/processed/summary_daily_z_WeStCOMS.rds"))



z_df <- map_dfr(c("WeStCOMS", "Linnhe", "Skye"), 
                ~readRDS(glue("{out_dir}/processed/summary_daily_z_{.x}.rds")) |>
                  mutate(region=.x))
