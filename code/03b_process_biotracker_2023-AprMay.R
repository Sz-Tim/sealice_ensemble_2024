# Compile output from long runs
# Sea lice ms 2024
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk

# These simulations are intended for hourly visualizations.
# Output files are created for each hour, with a column for each simulation and 
# a row for each occupied WeStCOMS element.



# setup -------------------------------------------------------------------

library(tidyverse); library(glue)
library(furrr)
library(carrier)
library(sf)
library(rstan)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
theme_set(theme_bw() + theme(panel.grid=element_blank()))

mesh_sf <- st_read("data/WeStCOMS2_mesh.gpkg") |> mutate(vol_top50m=area*pmin(depth, 50))
mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
mesh_i <- mesh_sf |> st_drop_geometry() |> select(i, area, vol_top50m)
out_dir <- "out/sim_2023-AprMay"
sim_i <- read_csv(glue("{out_dir}/sim_i.csv")) |>
  mutate(sim=paste0("sim_", i))
site_i <- read_csv("data/farm_sites_2023-AprMay.csv")
init_df <- full_join(
  site_i,
  read_csv("data/lice_daily_2023-AprMay.csv", skip=1,
             col_names=c("sepaSite", paste0("d_", 0:60)))
) |>
  pivot_longer(starts_with("d_"), names_to="date", values_to="density") |>
  mutate(date=ymd("2023-04-01") + as.numeric(str_sub(date, 3, -1))) |>
  left_join(read_csv("data/lice_biomass_2017-01-01_2023-12-31.csv") |>
              rename(fishTonnes=actualBiomassOnSiteTonnes))
dir.create(glue("{out_dir}/processed/hourly"), recursive=T, showWarnings=F)




# extract output ----------------------------------------------------------

# extract biotracker output: sim_[0-9][0-9].tar.gz
if(FALSE) {
  f_tgz <- dirf(out_dir, "tar.gz")
  walk(f_tgz, ~untar(.x, exdir=str_remove(.x, ".tar.gz")))
}




# ensemble weights --------------------------------------------------------

mod <- "pRE_huRE"
out_ensMixRE <- readRDS(glue("out/ensMix_3D_{mod}_stanfit.rds"))
dat_ensMixRE <- readRDS(glue("out/ensMix_3D_{mod}_standata.rds"))

iter_sub <- sample.int(length(rstan::extract(out_ensMixRE, pars="sigma")[[1]]), size=3000)
b_p_post <- rstan::extract(out_ensMixRE, pars="b_p")[[1]]


# particle densities ------------------------------------------------------

date_seq <- seq(ymd("2023-04-06"), ymd("2023-05-31"), by=1) |> str_remove_all("-")
ps_lims <- tibble(ens_mn=c(0,0),
                  ens_CL025=c(0,0),
                  ens_CL975=c(0,0),
                  ens_CI95width=c(0,0))
for(i in 1:length(date_seq)) {
  ps_i <- load_psteps_simSets(out_dir, mesh_i, sim_i, ncores=1, liceScale=1, 
                              stage=paste0("Mature_", date_seq[i]), per_m2=TRUE, trans="4th_rt")
  i_ts <- grep("^t_", names(ps_i), value=T)
  plan(multisession, workers=40)
  cat("Starting", date_seq[i])
  for(j in seq_along(i_ts)) {
    ps_j <- ps_i |> 
      filter(sim %in% dat_ensMixRE$sim_names) |>
      select("sim", "i", all_of(i_ts[j])) |>
      pivot_wider(names_from=sim, values_from=starts_with("t_")) 
    if(any(!is.na(ps_j$i))) {
      ps_j_mx <- ps_j |> select(all_of(dat_ensMixRE$sim_names)) |> as.matrix()
      ps_j_mx[is.na(ps_j_mx)] <- 0
      ps_j <- ps_j |>
        mutate(row=row_number()) |>
        mutate(ens=future_map(row, ~apply(b_p_post[iter_sub,], 1, function(x) sum(x * ps_j_mx[.x,]))),
               ens_mn=map_dbl(ens, mean),
               ens_CL025=map_dbl(ens, ~quantile(.x, probs=0.025)),
               ens_CL975=map_dbl(ens, ~quantile(.x, probs=0.975)),
               ens_CI95width=ens_CL975-ens_CL025) |>
        select(i, starts_with("ens_"))
      ps_j |>
        saveRDS(glue("{out_dir}/processed/hourly/Mature_{date_seq[i]}_{i_ts[j]}.rds"))
      ps_lims$ens_mn <- range(c(ps_lims$ens_mn, range(ps_j$ens_mn)))
      ps_lims$ens_CL025 <- range(c(ps_lims$ens_CL025, range(ps_j$ens_CL025)))
      ps_lims$ens_CL975 <- range(c(ps_lims$ens_CL975, range(ps_j$ens_CL975)))
      ps_lims$ens_CI95width <- range(c(ps_lims$ens_CI95width, range(ps_j$ens_CI95width))) 
    }
    cat("", j)
  }
  cat("\n")
  plan(sequential)
  gc()
}

# Averages: 4th root scale
ps_wide_rtrt <- readRDS(glue("{out_dir}/processed/psteps_wide_rtrt.rds"))

ps_long_rtrt <- ps_wide_rtrt |>
  pivot_longer(starts_with("t_"), names_to="date", values_to="rtrt_N_m2") |>
  mutate(date=first(init_df$date) + dhours(as.numeric(str_sub(date, 3, -1))))

tictoc::tic()
ps_avg_rtrt <- ps_long_rtrt |>
  calc_psteps_avg("rtrt_N_m2", ncores=10, mesh_sf=mesh_sf)
tictoc::toc()
saveRDS(ps_avg_rtrt, glue("{out_dir}/processed/psteps_avg_rtrt.rds"))

rm(list=ls() |> grep("ps_.*_rtrt", x=_, value=T))
gc()





# connectivity ------------------------------------------------------------

plan(multisession, workers=40)
c_long <- future_map_dfr(dirrf(out_dir, "connectivity.*csv"),
                         ~load_connectivity(.x, site_i$sepaSite, liceScale=1/24) |>
                           mutate(sim=str_sub(str_split_fixed(.x, "sim_", 2)[,2], 1, 2)))
plan(sequential)

# filter so that only sites with fish are included
active_df <- init_df |> filter(fishTonnes > 0) |> select(sepaSite, date)
c_long <- c_long |>
  inner_join(active_df, by=join_by(date, source==sepaSite)) |>
  inner_join(active_df, by=join_by(date, destination==sepaSite)) 

site_areas <- read_csv("data/farm_sites_100m_areas.csv") |> select(sepaSite, area)

# mean hourly IP for the week leading up to and including the date in column 'week'
c_weekly <- list(
  c_long |> 
    mutate(week=ceiling_date(date, "week")) |>
    calc_influx(destination, value, sim, week) |> 
    rename(sepaSite=destination) |>
    left_join(site_areas) |> 
    mutate(influx=influx/7,
           influx_m2=influx/area) |> 
    select(-area),
  c_long |> 
    mutate(week=ceiling_date(date, "week")) |>
    calc_outflux(source, value, sim, week, 
                 dest_areas=site_areas |> rename(destination=sepaSite)) |> 
    rename(sepaSite=source) |>
    mutate(outflux=outflux/7,
           outflux_m2=outflux_m2/7),
  c_long |> 
    mutate(week=ceiling_date(date, "week")) |>
    calc_self_infection(source, destination, value, sim, week) |> 
    rename(sepaSite=source) |>
    left_join(site_areas) |> 
    mutate(self=self/7,
           self_m2=self/area) |> 
    select(-area)
) |>
  reduce(full_join) |>
  left_join(init_df |> mutate(week=ceiling_date(date, "week")) |>
              group_by(sepaSite, week) |> summarise(fishTonnes=mean(fishTonnes, na.rm=T))) |>
  mutate(fishTonnes=if_else(is.na(fishTonnes), 0, fishTonnes)) |>
  filter(fishTonnes > 0) 

saveRDS(c_weekly, glue("{out_dir}/processed/connectivity_wk.rds"))

# mean hourly IP for each day
c_daily <- list(
  c_long |> 
    calc_influx(destination, value, sim, date) |> 
    rename(sepaSite=destination) |>
    left_join(site_areas) |>
    mutate(influx_m2=influx/area) |> 
    select(-area),
  c_long |> 
    calc_outflux(source, value, sim, date, 
                 dest_areas=site_areas |> rename(destination=sepaSite)) |> 
    rename(sepaSite=source),
  c_long |> 
    calc_self_infection(source, destination, value, sim, date) |> rename(sepaSite=source) |>
    left_join(site_areas) |> mutate(self_m2=self/area) |> select(-area)
) |>
  reduce(full_join) |>
  left_join(init_df) |>
  mutate(fishTonnes=if_else(is.na(fishTonnes), 0, fishTonnes)) |>
  filter(fishTonnes > 0)

saveRDS(c_daily, glue("{out_dir}/processed/connectivity_day.rds"))




