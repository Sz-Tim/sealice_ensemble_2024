# Long run simulations
# Sea lice ms 2024
# Tim Szewczyk
# tim.szewczyk@sams.ac.uk


# setup
library(tidyverse); library(glue)
library(sevcheck) # devtools::install_github("Sz-Tim/sevcheck")
library(biotrackR) # devtools::install_github("Sz-Tim/biotrackR")
library(doFuture)
theme_set(theme_bw() + theme(panel.grid=element_blank()))



# define parameters -------------------------------------------------------

cores_per_sim <- 12
parallel_sims <- 1
start_date <- "2023-03-17"
end_date <- "2023-05-31"
nDays <- length(seq(ymd(start_date), ymd(end_date), by=1))

os <- get_os()
dirs <- switch(
  get_os(),
  linux=list(proj=getwd(),
             mesh="/home/sa04ts/FVCOM_meshes",
             hydro="/media/archiver/common/sa01da-work/WeStCOMS2/Archive",
             jdk="/usr/local/java/jre1.8.0_211/bin/java",
             jar="/home/sa04ts/biotracker/biotracker.jar",
             out=glue("{getwd()}/out/sim_2023-MarMay")),
  windows=list(proj=getwd(),
               mesh="E:/hydroOut",
               hydro="E:/hydroOut/WeStCOMS2/Archive",
               jdk="C:/Users/sa04ts/.jdks/openjdk-22.0.1/bin/javaw",
               # jdk="C:/Users/sa04ts/.jdks/openjdk-19/bin/javaw",
               jar="C:/Users/sa04ts/OneDrive - SAMS/Projects/00_packages/biotracker/out/biotracker.jar",
               out="D:/sealice_ensembling/out/sim_2023-MarMay")
)
sim.i <- bind_rows(
  expand_grid(fixDepth="false",
              salinityThresh=c("high", "low"),
              salinityMort=c("false", "true"),
              eggTemp=c(F, T),
              swimSpeed=c(1e-4, 1e-3)),
  expand_grid(fixDepth="true",
              salinityThresh="high",
              salinityMort=c("false", "true"),
              eggTemp=c(F, T),
              swimSpeed=1e-4)
) |>
  mutate(i=str_pad(row_number(), 2, "left", "0"),
         outDir=glue("{dirs$out}/sim_{i}/"))
write_csv(sim.i, glue("{dirs$out}/sim_i.csv")) 
sim_seq <- 1:nrow(sim.i)



# set properties ----------------------------------------------------------

walk(sim_seq, ~dir.create(sim.i$outDir[.x], showWarnings=F))

walk(sim_seq, 
     ~set_biotracker_properties(
       properties_file_path=glue("{dirs$out}/sim_{sim.i$i[.x]}.properties"),
       parallelThreads=cores_per_sim,
       start_ymd=as.numeric(str_remove_all(start_date, "-")),
       # sitefile=glue("{dirs$proj}/data/farm_sites_2023-MarMay.csv"),
       # siteDensityPath=glue("{dirs$proj}/data/lice_daily_2023-MarMay.csv"),
       sitefile="D:/sealice_ensembling/data/farm_sites_2023-MarMay.csv",
       sitefileEnd="D:/sealice_ensembling/data/farm_sites_2023-MarMay.csv",
       siteDensityPath="D:/sealice_ensembling/data/lice_daily_2023-MarMay.csv",
       mesh1=glue("{dirs$mesh}/WeStCOMS2_mesh.nc"),
       datadir=glue("{dirs$hydro}/"),
       numberOfDays=nDays,
       nparts=10,
       checkOpenBoundaries="true",
       openBoundaryThresh=2000,
       fixDepth=sim.i$fixDepth[.x],
       variableDhV="false",
       salinityThreshMin=if_else(sim.i$salinityThresh[.x]=="high", 23, 20),
       salinityThreshMax=if_else(sim.i$salinityThresh[.x]=="high", 31, 23),
       salinityMort=sim.i$salinityMort[.x],
       eggTemp_b0=if_else(sim.i$eggTemp[.x], 0.17, 28.2),
       eggTemp_b1=if_else(sim.i$eggTemp[.x], 4.28, 0),
       vertSwimSpeedMean=-sim.i$swimSpeed[.x],
       vertSwimSpeedStd=sim.i$swimSpeed[.x]/5,
       sinkingRateMean=sim.i$swimSpeed[.x],
       sinkingRateStd=sim.i$swimSpeed[.x]/5,
       recordVertDistr="false",
       recordConnectivity="false",
       recordPsteps="true",
       splitPsteps="false",
       pstepsInterval=1,
       verboseSetUp="true"))


# run simulations ---------------------------------------------------------

# plan(multisession, workers=parallel_sims)
plan(sequential)
sim_seq <- 1:16
sim_sets <- split(sim_seq, rep(1:parallel_sims, length(sim_seq)/parallel_sims))
# foreach(j=1:parallel_sims) %dofuture% {
j <- 1
  for(i in sim_sets[[j]]) {
    setwd(dirs$proj)
    biotrackR::run_biotracker(
      jdk_path=dirs$jdk,
      jar_path=dirs$jar,
      f_properties=glue::glue("{dirs$out}/sim_{stringr::str_pad(i, 2, 'left', '0')}.properties"),
      sim_dir=glue::glue("{dirs$out}/sim_{stringr::str_pad(i, 2, 'left', '0')}/")
    )
  # }
}
plan(sequential)
