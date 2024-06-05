# Compile fish farm data
# Sea lice ms 2024
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


# setup
library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
library(zoo)
library(janitor)
library(sf)
library(raster)
library(spaths)
library(readxl)
conflicted::conflicts_prefer(dplyr::select, 
                             dplyr::filter, 
                             tidyr::extract)
theme_set(theme_bw() + theme(panel.grid=element_blank()))

dirf("code/fn") |> walk(source)



# load data ---------------------------------------------------------------

mesh.fp <- st_read("data/WeStCOMS2_meshFootprint.gpkg")

# Biomass is managed by SEPA, while sea lice is managed by MSS
# Separate site names are used, and there is not always total overlap
# SEPA sites are directly connected to licences, so they are used as the 
# locations, with the nearest sea lice counts applied to each

# Download data
download_lice_counts("2021-03-29", today(), "data/aquaculture_scot/ms_sea_lice_2021-present.csv")
download_fish_biomass("2017-01-01", today(), "data/aquaculture_scot/biomass_monthly_reports.csv")
download_ms_site_details("data/aquaculture_scot/ms_site_details.csv")
download_sepa_licenses("data/aquaculture_scot/se_licence_conditions.csv")

# Process datasets
sepa.locs <- read_csv("data/aquaculture_scot/se_licence_conditions.csv") |>
  filter(salmon==TRUE & waterType=="Seawater") |>
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) %>%
  filter(c(st_within(., mesh.fp, sparse=F))) |>
  rename(sepaSite=sepaSiteId) |>
  select(sepaSite, easting, northing, geometry)
mss.locs <- read_csv("data/aquaculture_scot/ms_site_details.csv") |>
  filter(aquacultureType=="Fish" & waterType=="Seawater") |>
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) %>%
  filter(c(st_within(., mesh.fp, sparse=F))) |>
  rename(mssID=marineScotlandSiteId) |>
  select(mssID, easting, northing, geometry)
sites.i <- sepa.locs %>%
  mutate(mssID=mss.locs$mssID[st_nearest_feature(., mss.locs)])

biomass_df <- read_csv("data/aquaculture_scot/biomass_monthly_reports.csv") |>
  mutate(year=year(date)) |>
  filter(waterType=="Seawater",
         between(date, ymd("2017-01-01"), ymd("2023-12-31"))) |>
  select(-siteName, -easting, -northing) |>
  inner_join(sites.i |> st_drop_geometry(), by="sepaSite") |>
  arrange(sepaSite, date)

# The 'CLEAN' xlsx has replaced <, >, etc:
# >X<Y = mean(x,y)
# <X = X/2
# <CoGP = 0.1
# >X = X+1
# 998 = 0
lice_df <- clean_mss_lice_xlsx("data/aquaculture_scot/Sea+Lice+data+2017+-+2020.xlsx") |>
  bind_rows(read_xlsx("data/aquaculture_scot/Sea+Lice+Publication+until+week+12%2C++Mar+2021_CLEAN.xlsx") |>
              janitor::clean_names(case="small_camel") |>
              rename(weeklyAverageAf=avgAfNoFish,
                     mitigation=mitigations) |>
              select(starts_with("site"), 
                     weekBeginning, weeklyAverageAf, mitigation) |>
              mutate(weekBeginning=date(weekBeginning))) |>
  bind_rows(read_csv("data/aquaculture_scot/ms_sea_lice_2021-present.csv") |>
              select(starts_with("site"), 
                     weekBeginning, weeklyAverageAf, mitigation)) |>
  filter(siteNo %in% sites.i$mssID) |>
  group_by(siteNo, weekBeginning) |>
  summarise(weeklyAverageAf=mean(weeklyAverageAf, na.rm=T),
            mitigation=paste(unique(mitigation), collapse=", ")) |>
  ungroup() |>
  mutate(mitigation=str_remove(mitigation, ", NA") |>
           str_remove(",$") |>
           str_replace(",,", ","),
         mitigation=na_if(mitigation, "NA")) |>
  right_join(sites.i |> st_drop_geometry() |> select(sepaSite, mssID), y=_, 
             by=c("mssID"="siteNo"), relationship="many-to-many")
write_csv(lice_df |> filter(!is.na(weeklyAverageAf)), "data/lice_data_nonInterpolated.csv")




# merge datasets ----------------------------------------------------------

combo.df <- full_join(lice_df |> 
                        mutate(year=year(weekBeginning), 
                               month=month(weekBeginning),
                               day=day(weekBeginning)) |>
                        select(sepaSite, 
                               year, month, day, 
                               weeklyAverageAf), 
                      biomass_df |> st_drop_geometry() |>
                        mutate(year=year(date),
                               month=month(date),
                               day=day(date)) |>
                        select(sepaSite,
                               year, month, day,
                               actualBiomassOnSiteTonnes) |>
                        full_join(expand_grid(sepaSite=unique(biomass_df$sepaSite),
                                              year=2017:2023,
                                              month=1:12,
                                              day=1)) |>
                        mutate(actualBiomassOnSiteTonnes=replace_na(actualBiomassOnSiteTonnes, 0))) |>
  mutate(week=week(ymd(paste(year, month, day, sep="-")))) |>
  group_by(sepaSite) |>
  filter(any(!is.na(actualBiomassOnSiteTonnes))) |>
  ungroup()

dates.df <- expand_grid(sepaSite=unique(combo.df$sepaSite),
                        date=seq(ymd("2017-01-01"), ymd("2023-12-31"), by=1)) |>
  mutate(year=year(date), month=month(date), day=day(date))

# Reporting was only required if averageAf exceeded a threshold:
#  - 2017.Jul - 2019.Jun.09: 3
#  - 2019.Jun.10 - 2021.Mar.28: 2
# https://www.gov.scot/publications/fish-health-inspectorate-sea-lice-information/
# 
# If no value was given, the smaller of median(day) and the median of (values 
# below the threshold from 2021-Apr to present) was used
default_3 <- median(filter(lice_df, weekBeginning > "2021-03-31" & weeklyAverageAf < 3)$weeklyAverageAf)
default_2 <- median(filter(lice_df, weekBeginning > "2021-03-31" & weeklyAverageAf < 2)$weeklyAverageAf)
# Assume 240 fish / tonne
out.df <- left_join(dates.df, combo.df, 
                    by=c("sepaSite", "year", "month", "day")) |>
  arrange(sepaSite, year, month, day) |>
  fill(week) |>
  # use reported lice value for the whole week
  group_by(sepaSite, year, week) |>
  mutate(weeklyAverageAf=median(weeklyAverageAf, na.rm=T)) |>
  # assume 0 biomass continues through the end of the month
  group_by(sepaSite, year, month) |>
  mutate(actualBiomassOnSiteTonnes=if_else(is.na(actualBiomassOnSiteTonnes) & 
                                             max(actualBiomassOnSiteTonnes, na.rm=T)==0, 
                                           0, 
                                           actualBiomassOnSiteTonnes)) |>
  # interpolate biomass to day (0 until first non-0 value, linear between values)
  group_by(sepaSite) |>
  mutate(actualBiomassOnSiteTonnes=as.numeric(na.fill(zoo(actualBiomassOnSiteTonnes), c(0, "extend", NA)))) |>
  # set biomass constant until end of final reported month, assume 0 after
  group_by(sepaSite, year, month) |>
  mutate(actualBiomassOnSiteTonnes=if_else(is.na(actualBiomassOnSiteTonnes), 
                                           median(actualBiomassOnSiteTonnes, na.rm=T),
                                           actualBiomassOnSiteTonnes)) |>
  ungroup() |>
  # if no lice counts, use min(global daily average, reporting threshold / 10)
  group_by(year, month, day) |>
  mutate(weeklyAverageAf=case_when(
    !is.na(weeklyAverageAf) ~ weeklyAverageAf,
    is.na(weeklyAverageAf) & date<"2019-06-10" ~ min(median(weeklyAverageAf, na.rm=T), default_3, na.rm=T),
    is.na(weeklyAverageAf) & date>"2019-06-09" ~ min(median(weeklyAverageAf, na.rm=T), default_2, na.rm=T))) |>
  ungroup() |>
  mutate(weeklyAverageAf=pmax(weeklyAverageAf, 0.01)) |>
  mutate(total_AF=actualBiomassOnSiteTonnes * 240 * weeklyAverageAf) |>
  arrange(sepaSite, date) |> 
  group_by(sepaSite) |>
  mutate(N_obs=sum(actualBiomassOnSiteTonnes > 0, na.rm=T)) |>
  filter(N_obs > 0) |>
  ungroup() |>
  filter(date <= rollforward(max(biomass_df$date)))


# save output -------------------------------------------------------------

out.df |>
  select(sepaSite, date, weeklyAverageAf, actualBiomassOnSiteTonnes, total_AF) |>
  write_csv(glue("data/lice_biomass_{min(out.df$date)}_{max(out.df$date)}.csv"))
out.df |>
  filter(date >= "2019-04-01") |>
  group_by(sepaSite) |>
  summarise() |>
  left_join(sites.i |> st_drop_geometry() |> select(sepaSite, easting, northing)) |>
  write_csv("data/farm_sites.csv")
out.df |>
  filter(date >= "2019-04-01") |>
  mutate(date.c=str_remove_all(date, "-")) |>
  select(sepaSite, date.c, total_AF) |>
  pivot_wider(names_from=date.c, values_from=total_AF) |>
  filter(sepaSite %in% read_csv("data/farm_sites.csv")$sepaSite) |>
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) |>
  write_csv(glue("data/lice_daily_2019-04-01_{max(out.df$date)}.csv"))

out.df |>
  filter(year(date) == 2023) |>
  group_by(sepaSite) |>
  summarise(active=any(actualBiomassOnSiteTonnes>0)) |>
  ungroup() |>
  filter(active) |>
  select(-active) |>
  left_join(sites.i |> st_drop_geometry() |> select(sepaSite, easting, northing)) |>
  write_csv(glue("data/farm_sites_2023.csv"))
out.df |>
  filter(year(date) == 2023) |>
  mutate(date.c=str_remove_all(date, "-")) |>
  select(sepaSite, date.c, total_AF) |>
  pivot_wider(names_from=date.c, values_from=total_AF) |>
  filter(sepaSite %in% read_csv("data/farm_sites_2023.csv")$sepaSite) |>
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) |>
  write_csv(glue("data/lice_daily_2023.csv"))

out.df |>
  filter(year(date) == 2023 & month(date) %in% c(4:5)) |>
  group_by(sepaSite) |>
  summarise(active=any(actualBiomassOnSiteTonnes>0)) |>
  ungroup() |>
  filter(active) |>
  select(-active) |>
  left_join(sites.i |> st_drop_geometry() |> select(sepaSite, easting, northing)) |>
  write_csv(glue("data/farm_sites_2023-AprMay.csv"))
out.df |>
  filter(year(date) == 2023 & month(date) %in% c(4:5)) |>
  mutate(date.c=str_remove_all(date, "-")) |>
  select(sepaSite, date.c, total_AF) |>
  pivot_wider(names_from=date.c, values_from=total_AF) |>
  filter(sepaSite %in% read_csv("data/farm_sites_2023-AprMay.csv")$sepaSite) |>
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) |>
  write_csv(glue("data/lice_daily_2023-AprMay.csv"))




# site characteristics ----------------------------------------------------

site_sf <- read_csv("data/farm_sites.csv") |> 
  st_as_sf(coords=c("easting", "northing"), crs=27700, remove=F) %>%
  mutate(fetch=raster::extract(raster::raster("data/log10_eu200m1a.tif"), ., 
                               small=T, buffer=1e2, fun=mean)) %>%
  mutate(fetch=if_else(is.na(fetch),
                       raster::extract(raster::raster("data/log10_eu200m1a.tif"), ., 
                                       small=T, buffer=5e2, fun=mean),
                       fetch))
write_csv(site_sf |> st_drop_geometry(), "data/farm_sites_fetch.csv")

# load ocean raster
mesh.r <- raster("data/ScotlandOcean_footprint.tif")
crs(mesh.r) <- CRS("+init=epsg:27700")
out_paths <- shortest_paths(mesh.r, 
                            site_sf)
write_csv(out_paths |> 
            as_tibble() |>
            mutate(origins=site_sf$sepaSite[origins],
                   destinations=site_sf$sepaSite[destinations]), 
          "data/site_pairwise_distances.csv")







# WeStCOMS WSPZs ----------------------------------------------------------

st_join(st_read("data/WeStCOMS2_mesh.gpkg") |> select(i, geom), 
        st_read("data/Wild_Salmonid_Protection_Zones_WGS84.gpkg") |> 
          st_transform(27700)) |>
  st_drop_geometry() |>
  write_csv("data/WeStCOMS2_WSPZ_lookup.csv")
