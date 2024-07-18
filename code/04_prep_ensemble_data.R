# Prepare ensemble data
# 
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk

# This script creates the dataset for fitting and validating the ensemble models.
# It reads the non-interpolated reported fish and lice data, identifies starts
# of production cycles, and subsets records from the start of each cycle to the 
# first application of a lice treatment. The daily infection pressure is then 
# read in for each simulation, followed by lagging, reduction due to daily
# mortality rates, and summed to estimate the resulting relative number of adult
# lice predicted by each simulation for each day. This is aligned with the data
# from the farms and stored for use in the next scripts.


# setup -------------------------------------------------------------------
library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck"); # dirrf() and get_lags()


site_i <- read_csv("data/farm_sites.csv") |>
  left_join(read_csv("data/farm_sites_100m_areas.csv") |> select(sepaSite, area))
init_df <- full_join(
  site_i,
  read_csv("data/lice_daily_2019-04-01_2023-12-31.csv", skip=1, 
           col_names=c("sepaSite", paste0("d_", 0:1736)))
) |>
  pivot_longer(starts_with("d_"), names_to="date", values_to="density") |>
  mutate(date=ymd("2019-04-01") + as.numeric(str_sub(date, 3, -1))) |>
  left_join(read_csv("data/lice_biomass_2017-01-01_2023-12-31.csv") |>
              rename(fishTonnes=actualBiomassOnSiteTonnes))
sim_i <- read_csv("out/sim_2019-2023/sim_i.csv") |>
  mutate(sim=paste0("sim_", i),
         lab_short=if_else(fixDepth, "2D", "3D")) |>
  group_by(lab_short) |>
  mutate(lab=paste0(lab_short, ".", row_number())) |>
  ungroup() |>
  select(sim, lab_short, lab) |>
  bind_rows(
    tibble(sim=c("predF1", "predF5", "predM", "predMRE", "sim_avg2D", "sim_avg3D", "null"),
           lab_short=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['M']", "Ens['Mix']", "Mean2D", "Mean3D", "Null"),
           lab=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['M']", "Ens['Mix']", "Mean2D", "Mean3D", "Null"))
  ) |>
  mutate(lab=factor(lab, 
                    levels=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['M']", "Ens['Mix']", "Mean2D", "Mean3D", 
                             paste0("2D.", 1:20), paste0("3D.", 1:20), "Null")),
         lab_short=factor(lab_short, 
                          levels=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['M']", "Ens['Mix']", "Mean2D", "Mean3D", 
                                   "2D", "3D", "Null")))
# sim_i <- dir("out/biotracker", "sim_i.csv", recursive=T, full.names=T) |>
#   map_dfr(~read_csv(.x) |> select(i, fixDepth) |> mutate(path=.x)) |>
#   filter((grepl("full", path) & i %in% paste0("0", 1:6)) |
#            (grepl("eggT", path) & i %in% str_pad(7:12, 2, "left", "0")) | 
#            (grepl("swim", path) & i %in% 13:20)) |>
#   arrange(i) |>
#   mutate(sim=paste0("sim_", i),
#          lab_short=if_else(fixDepth, "2D", "3D")) |>
#   group_by(lab_short) |>
#   mutate(lab=paste0(lab_short, ".", row_number())) |>
#   ungroup() |>
#   select(sim, lab_short, lab) |>
#   bind_rows(
#     tibble(sim=c("predF1", "predF5", "predM", "predMRE", "sim_avg2D", "sim_avg3D", "null"),
#            lab_short=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['M']", "Ens['Mix']", "Mean2D", "Mean3D", "Null"),
#            lab=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['M']", "Ens['Mix']", "Mean2D", "Mean3D", "Null"))
#   ) |>
#   mutate(lab=factor(lab, 
#                     levels=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['M']", "Ens['Mix']", "Mean2D", "Mean3D", 
#                              paste0("2D.", 1:20), paste0("3D.", 1:20), "Null")),
#          lab_short=factor(lab_short, 
#                           levels=c("Ens['Fc-1']", "Ens['Fc-5']", "Ens['M']", "Ens['Mix']", "Mean2D", "Mean3D", 
#                                    "2D", "3D", "Null")))


# load influx -------------------------------------------------------------

# read daily influx for all simulations
c_daily <- readRDS("out/sim_2019-2023/processed/connectivity_day.rds") |>
  mutate(sim=paste0("sim_", sim)) |>
  arrange(sim, sepaSite, date)
# c_daily <- dirrf("out/biotracker", "connectivity_day.rds") |>
#   map_dfr(~readRDS(.x) |> select(sepaSite, date, sim, influx_m2) |>
#             mutate(path=.x)) |>
#   rename(sim_og=sim) |>
#   # mutate(sim=paste0("sim_", str_split_fixed(path, "/", 4)[,3], "_", sim_og))
#   mutate(sim=paste0("sim_", sim_og)) |>
#   arrange(sim_og, sepaSite, date)



# identify production cycle starts to first lice treatments ---------------

newFarms_df <- init_df |>
  # for each farm, identify 'new' cycles by fishTonnes == 0 | fishTonnes < 25% of max
  group_by(sepaSite) |>
  mutate(fish_propMax=fishTonnes/max(fishTonnes, na.rm=T)) |>
  arrange(sepaSite, date) |>
  mutate(fishIn=!(is.na(fishTonnes) | fishTonnes==0 | fish_propMax < 0.25)) |>
  # find minimum of each period with 0-few fish
  group_by(sepaSite) |>
  mutate(fishInID=consecutive_id(fishIn)) |>
  filter(!fishIn) |>
  group_by(sepaSite, fishInID) |>
  slice_min(fishTonnes) |>
  slice_tail(n=1) |> # for consecutive 0s, select last day
  group_by(sepaSite) |>
  mutate(productionCycleNumber=row_number()) |>
  ungroup() |>
  select(sepaSite, date, productionCycleNumber) |>
  # join full dataset to fill in dates
  full_join(init_df, y=_) |>
  group_by(sepaSite) |>
  fill(productionCycleNumber) |>
  filter(!is.na(productionCycleNumber)) |>
  group_by(sepaSite, productionCycleNumber) |>
  ungroup() |>
  select(sepaSite, date, productionCycleNumber) |>
  # join with non-interpolated lice data
  inner_join(read_csv("data/lice_data_nonInterpolated.csv") |>
               rename(date=weekBeginning,
                      licePerFish=weeklyAverageAf) |>
               select(sepaSite, date, licePerFish, mitigation)) |>
  # identify first application of lice treatment (mitigation)
  group_by(sepaSite, productionCycleNumber) |>
  fill(mitigation, .direction="down") |>
  mutate(nDays=as.numeric(date - first(date))) |>
  ungroup() |>
  # take production cycle start until first treatment
  filter(is.na(mitigation))



# calculate relative expected adults from IP ------------------------------

maxLag <- 150 # presumed effective on-fish life span
surv <- c(1, rep(c(0.970, 0.998, 0.997, 0.942), times=c(15, 20, 20, maxLag-15-20-20)))
surv_ls <- map(1:length(surv), ~prod(surv[1:.x])) |>
  setNames(paste0("IP", 0:(length(surv)-1)))
valid_df <- c_daily |>
  select(sim, sepaSite, date, influx_m2) |>
  # fill in all dates for all sites (c_daily is sparse with 0's omitted)
  full_join(expand_grid(date=seq(ymd("2019-04-01"), ymd("2023-12-31"), by=1), 
                        sepaSite=unique(newFarms_df$sepaSite),
                        sim=unique(c_daily$sim))) |>
  mutate(influx_m2=replace_na(influx_m2, 0)) |>
  rename(IP=influx_m2) |>
  arrange(sim, sepaSite, date) |>
  # for each date, get lagged influx for up to 150 days 
  group_by(sim, sepaSite) |>
  get_lags(IP, n=maxLag) |>
  ungroup() |>
  # join only dates in newFarms_df
  inner_join(newFarms_df |> select(sepaSite, date, nDays)) |>
  rename(IP0=IP) |>
  # apply mortality for the corresponding number of days for each lag
  mutate(across(starts_with("IP"), ~.x*prod(surv_ls[[cur_column()]]))) |>
  # remove non-adult stages
  select(-(IP0:IP34)) |>
  pivot_longer(starts_with("IP")) |>
  # remove IP from before the start of the production cycle
  filter(as.numeric(str_remove(name, "IP")) <= nDays) |>
  # sum (mortality reduced) IP for all days = relative number of predicted adults
  group_by(sim, sepaSite, date) |>
  summarise(IP_adult=sum(value)) |>
  ungroup() |>
  inner_join(newFarms_df) |>
  filter(between(date, min(c_daily$date), max(c_daily$date)))



# final processing --------------------------------------------------------

valid_df <- valid_df |>
  # 4th root transformations due to right skew + 0's
  mutate(IP_adult_rtrt=IP_adult^0.25,
         licePerFish_rtrt=licePerFish^0.25) |>
  select(sim, sepaSite, productionCycleNumber, date, licePerFish_rtrt, IP_adult_rtrt) |>
  pivot_wider(names_from=sim, values_from=IP_adult_rtrt)|>
  # harmonic terms for day of year
  mutate(ydayCos=cos(yday(date)/366),
         ydaySin=sin(yday(date)/366),
         ydayCosXSin=ydayCos*ydaySin) |>
  # binary treatment by Code of Good Practice
  mutate(lice_g05=licePerFish_rtrt^4 > 0.5,
         year=year(date)) |>
  arrange(date, sepaSite) |>
  drop_na() |>
  mutate(rowNum=row_number(),
         sepaSiteNum=as.numeric(factor(sepaSite))) |>
  rowwise() |>
  mutate(sim_avg3D=mean(c_across(any_of(filter(sim_i, lab_short=="3D")$sim))),
         sim_avg2D=mean(c_across(any_of(filter(sim_i, lab_short=="2D")$sim)))) |>
  ungroup() |>
  mutate(across(starts_with("sim_"), ~.x - mean(.x), .names="c_{.col}"))

write_csv(valid_df, "out/valid_df.csv")
