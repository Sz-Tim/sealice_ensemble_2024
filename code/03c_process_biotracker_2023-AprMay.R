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

# mesh_sf <- st_read("data/WeStCOMS2_mesh.gpkg") |> mutate(vol_top50m=area*pmin(depth, 50))
# mesh_fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")
# mesh_i <- mesh_sf |> st_drop_geometry() |> select(i, area, vol_top50m)
out_dir <- "D:/sealice_ensembling/out/sim_2023-MarMay"
sim_i <- read_csv(glue("{out_dir}/sim_i.csv")) |>
  mutate(sim=paste0("sim_", i))
site_i <- read_csv("data/farm_sites_2023-MarMay.csv")
# init_df <- full_join(
#   site_i,
#   read_csv("data/lice_daily_2023-MarMay.csv", skip=1,
#              col_names=c("sepaSite", paste0("d_", 0:60)))
# ) |>
#   pivot_longer(starts_with("d_"), names_to="date", values_to="density") |>
#   mutate(date=ymd("2023-04-01") + as.numeric(str_sub(date, 3, -1))) |>
#   left_join(read_csv("data/lice_biomass_2017-01-01_2023-12-31.csv") |>
#               rename(fishTonnes=actualBiomassOnSiteTonnes))
# dir.create(glue("{out_dir}/processed/hourly"), recursive=T, showWarnings=F)




# extract output ----------------------------------------------------------

# extract biotracker output: sim_[0-9][0-9].tar.gz
if(FALSE) {
  f_tgz <- dirf(out_dir, "tar.gz")
  walk(f_tgz, ~untar(.x, exdir=str_remove(.x, ".tar.gz")))
}




# ensemble weights --------------------------------------------------------

mod <- "ranef"
out_ensMixRE <- readRDS(glue("out/ensembles/ensMix_3D_{mod}_stanfit.rds"))
dat_ensMixRE <- readRDS(glue("out/ensembles/ensMix_3D_{mod}_standata.rds"))

iter_sub <- sample.int(length(rstan::extract(out_ensMixRE, pars="sigma")[[1]]), size=3000)
b_p_post <- rstan::extract(out_ensMixRE, pars="b_p")[[1]]

rm(out_ensMixRE); gc()


# particle densities ------------------------------------------------------

date_seq <- seq(ymd("2023-05-03"), ymd("2023-05-31"), by=1) |> str_remove_all("-")
ps_lims <- tibble(ens_mn=c(0,0),
                  ens_sd=c(0,0),
                  ens_CL005=c(0,0),
                  ens_CL025=c(0,0),
                  ens_CL975=c(0,0),
                  ens_CL995=c(0,0),
                  ens_CI95width=c(0,0),
                  ens_CI99width=c(0,0))
for(i in 1:length(date_seq)) {
  ps_i <- load_psteps_simSets(out_dir, 
                              st_read("data/WeStCOMS2_mesh.gpkg") |>
                                mutate(vol_top50m=area*pmin(depth, 50)) |>
                                st_drop_geometry() |> 
                                select(i, area, vol_top50m), 
                              sim_i, ncores=1, liceScale=1, 
                              stage=paste0("Mature_", date_seq[i]), per_m2=TRUE, trans="4th_rt")
  i_ts <- grep("^t_", names(ps_i), value=T)
  plan(multisession, workers=6)
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
               ens_sd=map_dbl(ens, sd),
               ens_CL005=map_dbl(ens, ~quantile(.x, probs=0.005)),
               ens_CL025=map_dbl(ens, ~quantile(.x, probs=0.025)),
               ens_CL975=map_dbl(ens, ~quantile(.x, probs=0.975)),
               ens_CL995=map_dbl(ens, ~quantile(.x, probs=0.995)),
               ens_CI95width=ens_CL975-ens_CL025,
               ens_CI99width=ens_CL995-ens_CL005) |>
        select(i, starts_with("ens_"))
      gc()
      ps_j |>
        saveRDS(glue("{out_dir}/processed/hourly/Mature_{date_seq[i]}_{i_ts[j]}.rds"))
      ps_lims$ens_mn <- range(c(ps_lims$ens_mn, range(ps_j$ens_mn)))
      ps_lims$ens_sd <- range(c(ps_lims$ens_sd, range(ps_j$ens_sd)))
      ps_lims$ens_CL005 <- range(c(ps_lims$ens_CL005, range(ps_j$ens_CL005)))
      ps_lims$ens_CL025 <- range(c(ps_lims$ens_CL025, range(ps_j$ens_CL025)))
      ps_lims$ens_CL975 <- range(c(ps_lims$ens_CL975, range(ps_j$ens_CL975)))
      ps_lims$ens_CL995 <- range(c(ps_lims$ens_CL995, range(ps_j$ens_CL995)))
      ps_lims$ens_CI95width <- range(c(ps_lims$ens_CI95width, range(ps_j$ens_CI95width))) 
      ps_lims$ens_CI99width <- range(c(ps_lims$ens_CI99width, range(ps_j$ens_CI99width))) 
    }
    cat("", j)
  }
  cat("\n")
  plan(sequential)
  gc()
}
saveRDS(ps_lims, glue("{out_dir}/processed/hourly_Mature_pslims.rds"))





