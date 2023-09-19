# Title: Calculate center of gravity
# Authors: Rebecca Howard
# Date: 07/17/2023

## Libraries ----
library(dplyr)
library(here)
library(geosphere)
library(ggplot2)

## Functions ----
# Calculate center of gravity per time period
source(here('code/functions', 'COG_functions.R'))

## Pollock ----
### Eggs ----
# Load data
pkegg_hindcast <- readRDS(here('data', 'pk_egg_hindcast'))
df_pkegg1_cesm126 <- readRDS(here('data', 'df_pkegg1_cesm126.rds'))
df_pkegg2_cesm126 <- readRDS(here('data', 'df_pkegg2_cesm126.rds'))
df_pkegg3_cesm126 <- readRDS(here('data', 'df_pkegg3_cesm126.rds'))
df_pkegg1_cesm585 <- readRDS(here('data', 'df_pkegg1_cesm585.rds'))
df_pkegg2_cesm585 <- readRDS(here('data', 'df_pkegg2_cesm585.rds'))
df_pkegg3_cesm585 <- readRDS(here('data', 'df_pkegg3_cesm585.rds'))
df_pkegg1_gfdl126 <- readRDS(here('data', 'df_pkegg1_gfdl126.rds'))
df_pkegg2_gfdl126 <- readRDS(here('data', 'df_pkegg2_gfdl126.rds'))
df_pkegg3_gfdl126 <- readRDS(here('data', 'df_pkegg3_gfdl126.rds'))
df_pkegg1_gfdl585 <- readRDS(here('data', 'df_pkegg1_gfdl585.rds'))
df_pkegg2_gfdl585 <- readRDS(here('data', 'df_pkegg2_gfdl585.rds'))
df_pkegg3_gfdl585 <- readRDS(here('data', 'df_pkegg3_gfdl585.rds'))
df_pkegg1_miroc126 <- readRDS(here('data', 'df_pkegg1_miroc126.rds'))
df_pkegg2_miroc126 <- readRDS(here('data', 'df_pkegg2_miroc126.rds'))
df_pkegg3_miroc126 <- readRDS(here('data', 'df_pkegg3_miroc126.rds'))
df_pkegg1_miroc585 <- readRDS(here('data', 'df_pkegg1_miroc585.rds'))
df_pkegg2_miroc585 <- readRDS(here('data', 'df_pkegg2_miroc585.rds'))
df_pkegg3_miroc585 <- readRDS(here('data', 'df_pkegg3_miroc585.rds'))

pkegg1_sdm <- list(df_pkegg1_cesm126[[2]], df_pkegg1_cesm585[[2]],
                   df_pkegg1_gfdl126[[2]], df_pkegg1_gfdl585[[2]],
                   df_pkegg1_miroc126[[2]], df_pkegg1_miroc585[[2]])
pkegg1_temp <- list(df_pkegg1_cesm126[[3]], df_pkegg1_cesm585[[3]],
                    df_pkegg1_gfdl126[[3]], df_pkegg1_gfdl585[[3]],
                    df_pkegg1_miroc126[[3]], df_pkegg1_miroc585[[3]])

pkegg2_sdm <- list(df_pkegg2_cesm126[[2]], df_pkegg2_cesm585[[2]],
                   df_pkegg2_gfdl126[[2]], df_pkegg2_gfdl585[[2]],
                   df_pkegg2_miroc126[[2]], df_pkegg2_miroc585[[2]])
pkegg2_temp <- list(df_pkegg2_cesm126[[3]], df_pkegg2_cesm585[[3]],
                   df_pkegg2_gfdl126[[3]], df_pkegg2_gfdl585[[3]],
                   df_pkegg2_miroc126[[3]], df_pkegg2_miroc585[[3]])

pkegg3_sdm <- list(df_pkegg3_cesm126[[2]], df_pkegg3_cesm585[[2]],
                   df_pkegg3_gfdl126[[2]], df_pkegg3_gfdl585[[2]],
                   df_pkegg3_miroc126[[2]], df_pkegg3_miroc585[[2]])
pkegg3_temp <- list(df_pkegg3_cesm126[[3]], df_pkegg3_cesm585[[3]],
                   df_pkegg3_gfdl126[[3]], df_pkegg3_gfdl585[[3]],
                   df_pkegg3_miroc126[[3]], df_pkegg3_miroc585[[3]])

pkegg_COG <- COG_calc(pkegg_hindcast, pkegg1_sdm, pkegg2_sdm, pkegg3_sdm,
                      pkegg1_temp, pkegg2_temp, pkegg3_temp) 
saveRDS(pkegg_COG, here('data', 'pkegg_COG'))

rm(df_pkegg1_cesm126, df_pkegg1_cesm585, df_pkegg1_gfdl126, df_pkegg1_gfdl585,
   df_pkegg1_miroc126, df_pkegg1_miroc585, df_pkegg2_cesm126, df_pkegg2_cesm585, 
   df_pkegg2_gfdl126, df_pkegg2_gfdl585, df_pkegg2_miroc126, df_pkegg2_miroc585,
   df_pkegg3_cesm126, df_pkegg3_cesm585, df_pkegg3_gfdl126, df_pkegg3_gfdl585,
   df_pkegg3_miroc126, df_pkegg3_miroc585, pkegg1_sdm, pkegg1_temp, pkegg2_sdm,
   pkegg2_temp, pkegg3_sdm, pkegg3_temp)

### Larvae ----
# Load data
pklarvae_hindcast <- readRDS(here('data', 'pk_larvae_hindcast'))
df_pklarvae1_cesm126 <- readRDS(here('data', 'df_pklarvae1_cesm126.rds'))
df_pklarvae2_cesm126 <- readRDS(here('data', 'df_pklarvae2_cesm126.rds'))
df_pklarvae3_cesm126 <- readRDS(here('data', 'df_pklarvae3_cesm126.rds'))
df_pklarvae1_cesm585 <- readRDS(here('data', 'df_pklarvae1_cesm585.rds'))
df_pklarvae2_cesm585 <- readRDS(here('data', 'df_pklarvae2_cesm585.rds'))
df_pklarvae3_cesm585 <- readRDS(here('data', 'df_pklarvae3_cesm585.rds'))
df_pklarvae1_gfdl126 <- readRDS(here('data', 'df_pklarvae1_gfdl126.rds'))
df_pklarvae2_gfdl126 <- readRDS(here('data', 'df_pklarvae2_gfdl126.rds'))
df_pklarvae3_gfdl126 <- readRDS(here('data', 'df_pklarvae3_gfdl126.rds'))
df_pklarvae1_gfdl585 <- readRDS(here('data', 'df_pklarvae1_gfdl585.rds'))
df_pklarvae2_gfdl585 <- readRDS(here('data', 'df_pklarvae2_gfdl585.rds'))
df_pklarvae3_gfdl585 <- readRDS(here('data', 'df_pklarvae3_gfdl585.rds'))
df_pklarvae1_miroc126 <- readRDS(here('data', 'df_pklarvae1_miroc126.rds'))
df_pklarvae2_miroc126 <- readRDS(here('data', 'df_pklarvae2_miroc126.rds'))
df_pklarvae3_miroc126 <- readRDS(here('data', 'df_pklarvae3_miroc126.rds'))
df_pklarvae1_miroc585 <- readRDS(here('data', 'df_pklarvae1_miroc585.rds'))
df_pklarvae2_miroc585 <- readRDS(here('data', 'df_pklarvae2_miroc585.rds'))
df_pklarvae3_miroc585 <- readRDS(here('data', 'df_pklarvae3_miroc585.rds'))

pklarvae1_sdm <- list(df_pklarvae1_cesm126[[2]], df_pklarvae1_cesm585[[2]],
                      df_pklarvae1_gfdl126[[2]], df_pklarvae1_gfdl585[[2]],
                      df_pklarvae1_miroc126[[2]], df_pklarvae1_miroc585[[2]])
pklarvae1_temp <- list(df_pklarvae1_cesm126[[3]], df_pklarvae1_cesm585[[3]],
                       df_pklarvae1_gfdl126[[3]], df_pklarvae1_gfdl585[[3]],
                       df_pklarvae1_miroc126[[3]], df_pklarvae1_miroc585[[3]])

pklarvae2_sdm <- list(df_pklarvae2_cesm126[[2]], df_pklarvae2_cesm585[[2]],
                      df_pklarvae2_gfdl126[[2]], df_pklarvae2_gfdl585[[2]],
                      df_pklarvae2_miroc126[[2]], df_pklarvae2_miroc585[[2]])
pklarvae2_temp <- list(df_pklarvae2_cesm126[[3]], df_pklarvae2_cesm585[[3]],
                       df_pklarvae2_gfdl126[[3]], df_pklarvae2_gfdl585[[3]],
                       df_pklarvae2_miroc126[[3]], df_pklarvae2_miroc585[[3]])

pklarvae3_sdm <- list(df_pklarvae3_cesm126[[2]], df_pklarvae3_cesm585[[2]],
                      df_pklarvae3_gfdl126[[2]], df_pklarvae3_gfdl585[[2]],
                      df_pklarvae3_miroc126[[2]], df_pklarvae3_miroc585[[2]])
pklarvae3_temp <- list(df_pklarvae3_cesm126[[3]], df_pklarvae3_cesm585[[3]],
                       df_pklarvae3_gfdl126[[3]], df_pklarvae3_gfdl585[[3]],
                       df_pklarvae3_miroc126[[3]], df_pklarvae3_miroc585[[3]])

pklarvae_COG <- COG_calc(pklarvae_hindcast, pklarvae1_sdm, pklarvae2_sdm, pklarvae3_sdm,
                         pklarvae1_temp, pklarvae2_temp, pklarvae3_temp) 
saveRDS(pklarvae_COG, here('data', 'pklarvae_COG'))

rm(df_pklarvae1_cesm126, df_pklarvae1_cesm585, df_pklarvae1_gfdl126, df_pklarvae1_gfdl585,
   df_pklarvae1_miroc126, df_pklarvae1_miroc585, df_pklarvae2_cesm126, df_pklarvae2_cesm585, 
   df_pklarvae2_gfdl126, df_pklarvae2_gfdl585, df_pklarvae2_miroc126, df_pklarvae2_miroc585,
   df_pklarvae3_cesm126, df_pklarvae3_cesm585, df_pklarvae3_gfdl126, df_pklarvae3_gfdl585,
   df_pklarvae3_miroc126, df_pklarvae3_miroc585, pklarvae1_sdm, pklarvae1_temp, pklarvae2_sdm,
   pklarvae2_temp, pklarvae3_sdm, pklarvae3_temp)

## Flathead Sole ----
### Eggs ----
# Load data
fhsegg_hindcast <- readRDS(here('data', 'fhs_egg_hindcast'))
df_fhsegg1_cesm126 <- readRDS(here('data', 'df_fhsegg1_cesm126.rds'))
df_fhsegg2_cesm126 <- readRDS(here('data', 'df_fhsegg2_cesm126.rds'))
df_fhsegg3_cesm126 <- readRDS(here('data', 'df_fhsegg3_cesm126.rds'))
df_fhsegg1_cesm585 <- readRDS(here('data', 'df_fhsegg1_cesm585.rds'))
df_fhsegg2_cesm585 <- readRDS(here('data', 'df_fhsegg2_cesm585.rds'))
df_fhsegg3_cesm585 <- readRDS(here('data', 'df_fhsegg3_cesm585.rds'))
df_fhsegg1_gfdl126 <- readRDS(here('data', 'df_fhsegg1_gfdl126.rds'))
df_fhsegg2_gfdl126 <- readRDS(here('data', 'df_fhsegg2_gfdl126.rds'))
df_fhsegg3_gfdl126 <- readRDS(here('data', 'df_fhsegg3_gfdl126.rds'))
df_fhsegg1_gfdl585 <- readRDS(here('data', 'df_fhsegg1_gfdl585.rds'))
df_fhsegg2_gfdl585 <- readRDS(here('data', 'df_fhsegg2_gfdl585.rds'))
df_fhsegg3_gfdl585 <- readRDS(here('data', 'df_fhsegg3_gfdl585.rds'))
df_fhsegg1_miroc126 <- readRDS(here('data', 'df_fhsegg1_miroc126.rds'))
df_fhsegg2_miroc126 <- readRDS(here('data', 'df_fhsegg2_miroc126.rds'))
df_fhsegg3_miroc126 <- readRDS(here('data', 'df_fhsegg3_miroc126.rds'))
df_fhsegg1_miroc585 <- readRDS(here('data', 'df_fhsegg1_miroc585.rds'))
df_fhsegg2_miroc585 <- readRDS(here('data', 'df_fhsegg2_miroc585.rds'))
df_fhsegg3_miroc585 <- readRDS(here('data', 'df_fhsegg3_miroc585.rds'))

fhsegg1_sdm <- list(df_fhsegg1_cesm126[[2]], df_fhsegg1_cesm585[[2]],
                      df_fhsegg1_gfdl126[[2]], df_fhsegg1_gfdl585[[2]],
                      df_fhsegg1_miroc126[[2]], df_fhsegg1_miroc585[[2]])
fhsegg1_temp <- list(df_fhsegg1_cesm126[[3]], df_fhsegg1_cesm585[[3]],
                       df_fhsegg1_gfdl126[[3]], df_fhsegg1_gfdl585[[3]],
                       df_fhsegg1_miroc126[[3]], df_fhsegg1_miroc585[[3]])

fhsegg2_sdm <- list(df_fhsegg2_cesm126[[2]], df_fhsegg2_cesm585[[2]],
                      df_fhsegg2_gfdl126[[2]], df_fhsegg2_gfdl585[[2]],
                      df_fhsegg2_miroc126[[2]], df_fhsegg2_miroc585[[2]])
fhsegg2_temp <- list(df_fhsegg2_cesm126[[3]], df_fhsegg2_cesm585[[3]],
                       df_fhsegg2_gfdl126[[3]], df_fhsegg2_gfdl585[[3]],
                       df_fhsegg2_miroc126[[3]], df_fhsegg2_miroc585[[3]])

fhsegg3_sdm <- list(df_fhsegg3_cesm126[[2]], df_fhsegg3_cesm585[[2]],
                      df_fhsegg3_gfdl126[[2]], df_fhsegg3_gfdl585[[2]],
                      df_fhsegg3_miroc126[[2]], df_fhsegg3_miroc585[[2]])
fhsegg3_temp <- list(df_fhsegg3_cesm126[[3]], df_fhsegg3_cesm585[[3]],
                       df_fhsegg3_gfdl126[[3]], df_fhsegg3_gfdl585[[3]],
                       df_fhsegg3_miroc126[[3]], df_fhsegg3_miroc585[[3]])

fhsegg_COG <- COG_calc(fhsegg_hindcast, fhsegg1_sdm, fhsegg2_sdm, fhsegg3_sdm,
                         fhsegg1_temp, fhsegg2_temp, fhsegg3_temp) 
saveRDS(fhsegg_COG, here('data', 'fhsegg_COG'))

rm(df_fhsegg1_cesm126, df_fhsegg1_cesm585, df_fhsegg1_gfdl126, df_fhsegg1_gfdl585,
   df_fhsegg1_miroc126, df_fhsegg1_miroc585, df_fhsegg2_cesm126, df_fhsegg2_cesm585, 
   df_fhsegg2_gfdl126, df_fhsegg2_gfdl585, df_fhsegg2_miroc126, df_fhsegg2_miroc585,
   df_fhsegg3_cesm126, df_fhsegg3_cesm585, df_fhsegg3_gfdl126, df_fhsegg3_gfdl585,
   df_fhsegg3_miroc126, df_fhsegg3_miroc585, fhsegg1_sdm, fhsegg1_temp, fhsegg2_sdm,
   fhsegg2_temp, fhsegg3_sdm, fhsegg3_temp)

### Larvae ----
# Load data
fhslarvae_hindcast <- readRDS(here('data', 'fhs_larvae_hindcast'))
df_fhslarvae1_cesm126 <- readRDS(here('data', 'df_fhslarvae1_cesm126.rds'))
df_fhslarvae2_cesm126 <- readRDS(here('data', 'df_fhslarvae2_cesm126.rds'))
df_fhslarvae3_cesm126 <- readRDS(here('data', 'df_fhslarvae3_cesm126.rds'))
df_fhslarvae1_cesm585 <- readRDS(here('data', 'df_fhslarvae1_cesm585.rds'))
df_fhslarvae2_cesm585 <- readRDS(here('data', 'df_fhslarvae2_cesm585.rds'))
df_fhslarvae3_cesm585 <- readRDS(here('data', 'df_fhslarvae3_cesm585.rds'))
df_fhslarvae1_gfdl126 <- readRDS(here('data', 'df_fhslarvae1_gfdl126.rds'))
df_fhslarvae2_gfdl126 <- readRDS(here('data', 'df_fhslarvae2_gfdl126.rds'))
df_fhslarvae3_gfdl126 <- readRDS(here('data', 'df_fhslarvae3_gfdl126.rds'))
df_fhslarvae1_gfdl585 <- readRDS(here('data', 'df_fhslarvae1_gfdl585.rds'))
df_fhslarvae2_gfdl585 <- readRDS(here('data', 'df_fhslarvae2_gfdl585.rds'))
df_fhslarvae3_gfdl585 <- readRDS(here('data', 'df_fhslarvae3_gfdl585.rds'))
df_fhslarvae1_miroc126 <- readRDS(here('data', 'df_fhslarvae1_miroc126.rds'))
df_fhslarvae2_miroc126 <- readRDS(here('data', 'df_fhslarvae2_miroc126.rds'))
df_fhslarvae3_miroc126 <- readRDS(here('data', 'df_fhslarvae3_miroc126.rds'))
df_fhslarvae1_miroc585 <- readRDS(here('data', 'df_fhslarvae1_miroc585.rds'))
df_fhslarvae2_miroc585 <- readRDS(here('data', 'df_fhslarvae2_miroc585.rds'))
df_fhslarvae3_miroc585 <- readRDS(here('data', 'df_fhslarvae3_miroc585.rds'))

fhslarvae1_sdm <- list(df_fhslarvae1_cesm126[[2]], df_fhslarvae1_cesm585[[2]],
                    df_fhslarvae1_gfdl126[[2]], df_fhslarvae1_gfdl585[[2]],
                    df_fhslarvae1_miroc126[[2]], df_fhslarvae1_miroc585[[2]])
fhslarvae1_temp <- list(df_fhslarvae1_cesm126[[3]], df_fhslarvae1_cesm585[[3]],
                     df_fhslarvae1_gfdl126[[3]], df_fhslarvae1_gfdl585[[3]],
                     df_fhslarvae1_miroc126[[3]], df_fhslarvae1_miroc585[[3]])

fhslarvae2_sdm <- list(df_fhslarvae2_cesm126[[2]], df_fhslarvae2_cesm585[[2]],
                    df_fhslarvae2_gfdl126[[2]], df_fhslarvae2_gfdl585[[2]],
                    df_fhslarvae2_miroc126[[2]], df_fhslarvae2_miroc585[[2]])
fhslarvae2_temp <- list(df_fhslarvae2_cesm126[[3]], df_fhslarvae2_cesm585[[3]],
                     df_fhslarvae2_gfdl126[[3]], df_fhslarvae2_gfdl585[[3]],
                     df_fhslarvae2_miroc126[[3]], df_fhslarvae2_miroc585[[3]])

fhslarvae3_sdm <- list(df_fhslarvae3_cesm126[[2]], df_fhslarvae3_cesm585[[2]],
                    df_fhslarvae3_gfdl126[[2]], df_fhslarvae3_gfdl585[[2]],
                    df_fhslarvae3_miroc126[[2]], df_fhslarvae3_miroc585[[2]])
fhslarvae3_temp <- list(df_fhslarvae3_cesm126[[3]], df_fhslarvae3_cesm585[[3]],
                     df_fhslarvae3_gfdl126[[3]], df_fhslarvae3_gfdl585[[3]],
                     df_fhslarvae3_miroc126[[3]], df_fhslarvae3_miroc585[[3]])

fhslarvae_COG <- COG_calc(fhslarvae_hindcast, fhslarvae1_sdm, fhslarvae2_sdm, fhslarvae3_sdm,
                       fhslarvae1_temp, fhslarvae2_temp, fhslarvae3_temp) 
saveRDS(fhslarvae_COG, here('data', 'fhslarvae_COG'))

rm(df_fhslarvae1_cesm126, df_fhslarvae1_cesm585, df_fhslarvae1_gfdl126, df_fhslarvae1_gfdl585,
   df_fhslarvae1_miroc126, df_fhslarvae1_miroc585, df_fhslarvae2_cesm126, df_fhslarvae2_cesm585, 
   df_fhslarvae2_gfdl126, df_fhslarvae2_gfdl585, df_fhslarvae2_miroc126, df_fhslarvae2_miroc585,
   df_fhslarvae3_cesm126, df_fhslarvae3_cesm585, df_fhslarvae3_gfdl126, df_fhslarvae3_gfdl585,
   df_fhslarvae3_miroc126, df_fhslarvae3_miroc585, fhslarvae1_sdm, fhslarvae1_temp, fhslarvae2_sdm,
   fhslarvae2_temp, fhslarvae3_sdm, fhslarvae3_temp)

## Alaska Plaice ----
### Eggs ----
# Load data
akpegg_hindcast <- readRDS(here('data', 'akp_egg_hindcast'))
df_akpegg1_cesm126 <- readRDS(here('data', 'df_akpegg1_cesm126.rds'))
df_akpegg2_cesm126 <- readRDS(here('data', 'df_akpegg2_cesm126.rds'))
df_akpegg3_cesm126 <- readRDS(here('data', 'df_akpegg3_cesm126.rds'))
df_akpegg1_cesm585 <- readRDS(here('data', 'df_akpegg1_cesm585.rds'))
df_akpegg2_cesm585 <- readRDS(here('data', 'df_akpegg2_cesm585.rds'))
df_akpegg3_cesm585 <- readRDS(here('data', 'df_akpegg3_cesm585.rds'))
df_akpegg1_gfdl126 <- readRDS(here('data', 'df_akpegg1_gfdl126.rds'))
df_akpegg2_gfdl126 <- readRDS(here('data', 'df_akpegg2_gfdl126.rds'))
df_akpegg3_gfdl126 <- readRDS(here('data', 'df_akpegg3_gfdl126.rds'))
df_akpegg1_gfdl585 <- readRDS(here('data', 'df_akpegg1_gfdl585.rds'))
df_akpegg2_gfdl585 <- readRDS(here('data', 'df_akpegg2_gfdl585.rds'))
df_akpegg3_gfdl585 <- readRDS(here('data', 'df_akpegg3_gfdl585.rds'))
df_akpegg1_miroc126 <- readRDS(here('data', 'df_akpegg1_miroc126.rds'))
df_akpegg2_miroc126 <- readRDS(here('data', 'df_akpegg2_miroc126.rds'))
df_akpegg3_miroc126 <- readRDS(here('data', 'df_akpegg3_miroc126.rds'))
df_akpegg1_miroc585 <- readRDS(here('data', 'df_akpegg1_miroc585.rds'))
df_akpegg2_miroc585 <- readRDS(here('data', 'df_akpegg2_miroc585.rds'))
df_akpegg3_miroc585 <- readRDS(here('data', 'df_akpegg3_miroc585.rds'))

akpegg1_sdm <- list(df_akpegg1_cesm126[[2]], df_akpegg1_cesm585[[2]],
                    df_akpegg1_gfdl126[[2]], df_akpegg1_gfdl585[[2]],
                    df_akpegg1_miroc126[[2]], df_akpegg1_miroc585[[2]])
akpegg1_temp <- list(df_akpegg1_cesm126[[3]], df_akpegg1_cesm585[[3]],
                     df_akpegg1_gfdl126[[3]], df_akpegg1_gfdl585[[3]],
                     df_akpegg1_miroc126[[3]], df_akpegg1_miroc585[[3]])

akpegg2_sdm <- list(df_akpegg2_cesm126[[2]], df_akpegg2_cesm585[[2]],
                    df_akpegg2_gfdl126[[2]], df_akpegg2_gfdl585[[2]],
                    df_akpegg2_miroc126[[2]], df_akpegg2_miroc585[[2]])
akpegg2_temp <- list(df_akpegg2_cesm126[[3]], df_akpegg2_cesm585[[3]],
                     df_akpegg2_gfdl126[[3]], df_akpegg2_gfdl585[[3]],
                     df_akpegg2_miroc126[[3]], df_akpegg2_miroc585[[3]])

akpegg3_sdm <- list(df_akpegg3_cesm126[[2]], df_akpegg3_cesm585[[2]],
                    df_akpegg3_gfdl126[[2]], df_akpegg3_gfdl585[[2]],
                    df_akpegg3_miroc126[[2]], df_akpegg3_miroc585[[2]])
akpegg3_temp <- list(df_akpegg3_cesm126[[3]], df_akpegg3_cesm585[[3]],
                     df_akpegg3_gfdl126[[3]], df_akpegg3_gfdl585[[3]],
                     df_akpegg3_miroc126[[3]], df_akpegg3_miroc585[[3]])

akpegg_COG <- COG_calc(akpegg_hindcast, akpegg1_sdm, akpegg2_sdm, akpegg3_sdm,
                       akpegg1_temp, akpegg2_temp, akpegg3_temp) 
saveRDS(akpegg_COG, here('data', 'akpegg_COG'))

rm(df_akpegg1_cesm126, df_akpegg1_cesm585, df_akpegg1_gfdl126, df_akpegg1_gfdl585,
   df_akpegg1_miroc126, df_akpegg1_miroc585, df_akpegg2_cesm126, df_akpegg2_cesm585, 
   df_akpegg2_gfdl126, df_akpegg2_gfdl585, df_akpegg2_miroc126, df_akpegg2_miroc585,
   df_akpegg3_cesm126, df_akpegg3_cesm585, df_akpegg3_gfdl126, df_akpegg3_gfdl585,
   df_akpegg3_miroc126, df_akpegg3_miroc585, akpegg1_sdm, akpegg1_temp, akpegg2_sdm,
   akpegg2_temp, akpegg3_sdm, akpegg3_temp)

### Larvae ----
# Load data
akplarvae_hindcast <- readRDS(here('data', 'akp_larvae_hindcast'))
df_akplarvae1_cesm126 <- readRDS(here('data', 'df_akplarvae1_cesm126.rds'))
df_akplarvae2_cesm126 <- readRDS(here('data', 'df_akplarvae2_cesm126.rds'))
df_akplarvae3_cesm126 <- readRDS(here('data', 'df_akplarvae3_cesm126.rds'))
df_akplarvae1_cesm585 <- readRDS(here('data', 'df_akplarvae1_cesm585.rds'))
df_akplarvae2_cesm585 <- readRDS(here('data', 'df_akplarvae2_cesm585.rds'))
df_akplarvae3_cesm585 <- readRDS(here('data', 'df_akplarvae3_cesm585.rds'))
df_akplarvae1_gfdl126 <- readRDS(here('data', 'df_akplarvae1_gfdl126.rds'))
df_akplarvae2_gfdl126 <- readRDS(here('data', 'df_akplarvae2_gfdl126.rds'))
df_akplarvae3_gfdl126 <- readRDS(here('data', 'df_akplarvae3_gfdl126.rds'))
df_akplarvae1_gfdl585 <- readRDS(here('data', 'df_akplarvae1_gfdl585.rds'))
df_akplarvae2_gfdl585 <- readRDS(here('data', 'df_akplarvae2_gfdl585.rds'))
df_akplarvae3_gfdl585 <- readRDS(here('data', 'df_akplarvae3_gfdl585.rds'))
df_akplarvae1_miroc126 <- readRDS(here('data', 'df_akplarvae1_miroc126.rds'))
df_akplarvae2_miroc126 <- readRDS(here('data', 'df_akplarvae2_miroc126.rds'))
df_akplarvae3_miroc126 <- readRDS(here('data', 'df_akplarvae3_miroc126.rds'))
df_akplarvae1_miroc585 <- readRDS(here('data', 'df_akplarvae1_miroc585.rds'))
df_akplarvae2_miroc585 <- readRDS(here('data', 'df_akplarvae2_miroc585.rds'))
df_akplarvae3_miroc585 <- readRDS(here('data', 'df_akplarvae3_miroc585.rds'))

akplarvae1_sdm <- list(df_akplarvae1_cesm126[[2]], df_akplarvae1_cesm585[[2]],
                    df_akplarvae1_gfdl126[[2]], df_akplarvae1_gfdl585[[2]],
                    df_akplarvae1_miroc126[[2]], df_akplarvae1_miroc585[[2]])
akplarvae1_temp <- list(df_akplarvae1_cesm126[[3]], df_akplarvae1_cesm585[[3]],
                     df_akplarvae1_gfdl126[[3]], df_akplarvae1_gfdl585[[3]],
                     df_akplarvae1_miroc126[[3]], df_akplarvae1_miroc585[[3]])

akplarvae2_sdm <- list(df_akplarvae2_cesm126[[2]], df_akplarvae2_cesm585[[2]],
                    df_akplarvae2_gfdl126[[2]], df_akplarvae2_gfdl585[[2]],
                    df_akplarvae2_miroc126[[2]], df_akplarvae2_miroc585[[2]])
akplarvae2_temp <- list(df_akplarvae2_cesm126[[3]], df_akplarvae2_cesm585[[3]],
                     df_akplarvae2_gfdl126[[3]], df_akplarvae2_gfdl585[[3]],
                     df_akplarvae2_miroc126[[3]], df_akplarvae2_miroc585[[3]])

akplarvae3_sdm <- list(df_akplarvae3_cesm126[[2]], df_akplarvae3_cesm585[[2]],
                    df_akplarvae3_gfdl126[[2]], df_akplarvae3_gfdl585[[2]],
                    df_akplarvae3_miroc126[[2]], df_akplarvae3_miroc585[[2]])
akplarvae3_temp <- list(df_akplarvae3_cesm126[[3]], df_akplarvae3_cesm585[[3]],
                     df_akplarvae3_gfdl126[[3]], df_akplarvae3_gfdl585[[3]],
                     df_akplarvae3_miroc126[[3]], df_akplarvae3_miroc585[[3]])

akplarvae_COG <- COG_calc(akplarvae_hindcast, akplarvae1_sdm, akplarvae2_sdm, akplarvae3_sdm,
                       akplarvae1_temp, akplarvae2_temp, akplarvae3_temp) 
saveRDS(akplarvae_COG, here('data', 'akplarvae_COG'))

rm(df_akplarvae1_cesm126, df_akplarvae1_cesm585, df_akplarvae1_gfdl126, df_akplarvae1_gfdl585,
   df_akplarvae1_miroc126, df_akplarvae1_miroc585, df_akplarvae2_cesm126, df_akplarvae2_cesm585, 
   df_akplarvae2_gfdl126, df_akplarvae2_gfdl585, df_akplarvae2_miroc126, df_akplarvae2_miroc585,
   df_akplarvae3_cesm126, df_akplarvae3_cesm585, df_akplarvae3_gfdl126, df_akplarvae3_gfdl585,
   df_akplarvae3_miroc126, df_akplarvae3_miroc585, akplarvae1_sdm, akplarvae1_temp, akplarvae2_sdm,
   akplarvae2_temp, akplarvae3_sdm, akplarvae3_temp)

## Yellowfin Sole ----
### Larvae ----
# Load data
yfslarvae_hindcast <- readRDS(here('data', 'yfs_larvae_hindcast'))
df_yfslarvae1_cesm126 <- readRDS(here('data', 'df_yfslarvae1_cesm126.rds'))
df_yfslarvae2_cesm126 <- readRDS(here('data', 'df_yfslarvae2_cesm126.rds'))
df_yfslarvae3_cesm126 <- readRDS(here('data', 'df_yfslarvae3_cesm126.rds'))
df_yfslarvae1_cesm585 <- readRDS(here('data', 'df_yfslarvae1_cesm585.rds'))
df_yfslarvae2_cesm585 <- readRDS(here('data', 'df_yfslarvae2_cesm585.rds'))
df_yfslarvae3_cesm585 <- readRDS(here('data', 'df_yfslarvae3_cesm585.rds'))
df_yfslarvae1_gfdl126 <- readRDS(here('data', 'df_yfslarvae1_gfdl126.rds'))
df_yfslarvae2_gfdl126 <- readRDS(here('data', 'df_yfslarvae2_gfdl126.rds'))
df_yfslarvae3_gfdl126 <- readRDS(here('data', 'df_yfslarvae3_gfdl126.rds'))
df_yfslarvae1_gfdl585 <- readRDS(here('data', 'df_yfslarvae1_gfdl585.rds'))
df_yfslarvae2_gfdl585 <- readRDS(here('data', 'df_yfslarvae2_gfdl585.rds'))
df_yfslarvae3_gfdl585 <- readRDS(here('data', 'df_yfslarvae3_gfdl585.rds'))
df_yfslarvae1_miroc126 <- readRDS(here('data', 'df_yfslarvae1_miroc126.rds'))
df_yfslarvae2_miroc126 <- readRDS(here('data', 'df_yfslarvae2_miroc126.rds'))
df_yfslarvae3_miroc126 <- readRDS(here('data', 'df_yfslarvae3_miroc126.rds'))
df_yfslarvae1_miroc585 <- readRDS(here('data', 'df_yfslarvae1_miroc585.rds'))
df_yfslarvae2_miroc585 <- readRDS(here('data', 'df_yfslarvae2_miroc585.rds'))
df_yfslarvae3_miroc585 <- readRDS(here('data', 'df_yfslarvae3_miroc585.rds'))

yfslarvae1_sdm <- list(df_yfslarvae1_cesm126[[2]], df_yfslarvae1_cesm585[[2]],
                       df_yfslarvae1_gfdl126[[2]], df_yfslarvae1_gfdl585[[2]],
                       df_yfslarvae1_miroc126[[2]], df_yfslarvae1_miroc585[[2]])
yfslarvae1_temp <- list(df_yfslarvae1_cesm126[[3]], df_yfslarvae1_cesm585[[3]],
                        df_yfslarvae1_gfdl126[[3]], df_yfslarvae1_gfdl585[[3]],
                        df_yfslarvae1_miroc126[[3]], df_yfslarvae1_miroc585[[3]])

yfslarvae2_sdm <- list(df_yfslarvae2_cesm126[[2]], df_yfslarvae2_cesm585[[2]],
                       df_yfslarvae2_gfdl126[[2]], df_yfslarvae2_gfdl585[[2]],
                       df_yfslarvae2_miroc126[[2]], df_yfslarvae2_miroc585[[2]])
yfslarvae2_temp <- list(df_yfslarvae2_cesm126[[3]], df_yfslarvae2_cesm585[[3]],
                        df_yfslarvae2_gfdl126[[3]], df_yfslarvae2_gfdl585[[3]],
                        df_yfslarvae2_miroc126[[3]], df_yfslarvae2_miroc585[[3]])

yfslarvae3_sdm <- list(df_yfslarvae3_cesm126[[2]], df_yfslarvae3_cesm585[[2]],
                       df_yfslarvae3_gfdl126[[2]], df_yfslarvae3_gfdl585[[2]],
                       df_yfslarvae3_miroc126[[2]], df_yfslarvae3_miroc585[[2]])
yfslarvae3_temp <- list(df_yfslarvae3_cesm126[[3]], df_yfslarvae3_cesm585[[3]],
                        df_yfslarvae3_gfdl126[[3]], df_yfslarvae3_gfdl585[[3]],
                        df_yfslarvae3_miroc126[[3]], df_yfslarvae3_miroc585[[3]])

yfslarvae_COG <- COG_calc(yfslarvae_hindcast, yfslarvae1_sdm, yfslarvae2_sdm, yfslarvae3_sdm,
                          yfslarvae1_temp, yfslarvae2_temp, yfslarvae3_temp) 
saveRDS(yfslarvae_COG, here('data', 'yfslarvae_COG'))

rm(df_yfslarvae1_cesm126, df_yfslarvae1_cesm585, df_yfslarvae1_gfdl126, df_yfslarvae1_gfdl585,
   df_yfslarvae1_miroc126, df_yfslarvae1_miroc585, df_yfslarvae2_cesm126, df_yfslarvae2_cesm585, 
   df_yfslarvae2_gfdl126, df_yfslarvae2_gfdl585, df_yfslarvae2_miroc126, df_yfslarvae2_miroc585,
   df_yfslarvae3_cesm126, df_yfslarvae3_cesm585, df_yfslarvae3_gfdl126, df_yfslarvae3_gfdl585,
   df_yfslarvae3_miroc126, df_yfslarvae3_miroc585, yfslarvae1_sdm, yfslarvae1_temp, yfslarvae2_sdm,
   yfslarvae2_temp, yfslarvae3_sdm, yfslarvae3_temp)

## Northern Rock Sole ----
### Larvae ----
# Load data
nrslarvae_hindcast <- readRDS(here('data', 'nrs_larvae_hindcast'))
df_nrslarvae1_cesm126 <- readRDS(here('data', 'df_nrslarvae1_cesm126.rds'))
df_nrslarvae2_cesm126 <- readRDS(here('data', 'df_nrslarvae2_cesm126.rds'))
df_nrslarvae3_cesm126 <- readRDS(here('data', 'df_nrslarvae3_cesm126.rds'))
df_nrslarvae1_cesm585 <- readRDS(here('data', 'df_nrslarvae1_cesm585.rds'))
df_nrslarvae2_cesm585 <- readRDS(here('data', 'df_nrslarvae2_cesm585.rds'))
df_nrslarvae3_cesm585 <- readRDS(here('data', 'df_nrslarvae3_cesm585.rds'))
df_nrslarvae1_gfdl126 <- readRDS(here('data', 'df_nrslarvae1_gfdl126.rds'))
df_nrslarvae2_gfdl126 <- readRDS(here('data', 'df_nrslarvae2_gfdl126.rds'))
df_nrslarvae3_gfdl126 <- readRDS(here('data', 'df_nrslarvae3_gfdl126.rds'))
df_nrslarvae1_gfdl585 <- readRDS(here('data', 'df_nrslarvae1_gfdl585.rds'))
df_nrslarvae2_gfdl585 <- readRDS(here('data', 'df_nrslarvae2_gfdl585.rds'))
df_nrslarvae3_gfdl585 <- readRDS(here('data', 'df_nrslarvae3_gfdl585.rds'))
df_nrslarvae1_miroc126 <- readRDS(here('data', 'df_nrslarvae1_miroc126.rds'))
df_nrslarvae2_miroc126 <- readRDS(here('data', 'df_nrslarvae2_miroc126.rds'))
df_nrslarvae3_miroc126 <- readRDS(here('data', 'df_nrslarvae3_miroc126.rds'))
df_nrslarvae1_miroc585 <- readRDS(here('data', 'df_nrslarvae1_miroc585.rds'))
df_nrslarvae2_miroc585 <- readRDS(here('data', 'df_nrslarvae2_miroc585.rds'))
df_nrslarvae3_miroc585 <- readRDS(here('data', 'df_nrslarvae3_miroc585.rds'))

nrslarvae1_sdm <- list(df_nrslarvae1_cesm126[[2]], df_nrslarvae1_cesm585[[2]],
                       df_nrslarvae1_gfdl126[[2]], df_nrslarvae1_gfdl585[[2]],
                       df_nrslarvae1_miroc126[[2]], df_nrslarvae1_miroc585[[2]])
nrslarvae1_temp <- list(df_nrslarvae1_cesm126[[3]], df_nrslarvae1_cesm585[[3]],
                        df_nrslarvae1_gfdl126[[3]], df_nrslarvae1_gfdl585[[3]],
                        df_nrslarvae1_miroc126[[3]], df_nrslarvae1_miroc585[[3]])

nrslarvae2_sdm <- list(df_nrslarvae2_cesm126[[2]], df_nrslarvae2_cesm585[[2]],
                       df_nrslarvae2_gfdl126[[2]], df_nrslarvae2_gfdl585[[2]],
                       df_nrslarvae2_miroc126[[2]], df_nrslarvae2_miroc585[[2]])
nrslarvae2_temp <- list(df_nrslarvae2_cesm126[[3]], df_nrslarvae2_cesm585[[3]],
                        df_nrslarvae2_gfdl126[[3]], df_nrslarvae2_gfdl585[[3]],
                        df_nrslarvae2_miroc126[[3]], df_nrslarvae2_miroc585[[3]])

nrslarvae3_sdm <- list(df_nrslarvae3_cesm126[[2]], df_nrslarvae3_cesm585[[2]],
                       df_nrslarvae3_gfdl126[[2]], df_nrslarvae3_gfdl585[[2]],
                       df_nrslarvae3_miroc126[[2]], df_nrslarvae3_miroc585[[2]])
nrslarvae3_temp <- list(df_nrslarvae3_cesm126[[3]], df_nrslarvae3_cesm585[[3]],
                        df_nrslarvae3_gfdl126[[3]], df_nrslarvae3_gfdl585[[3]],
                        df_nrslarvae3_miroc126[[3]], df_nrslarvae3_miroc585[[3]])

nrslarvae_COG <- COG_calc(nrslarvae_hindcast, nrslarvae1_sdm, nrslarvae2_sdm, nrslarvae3_sdm,
                          nrslarvae1_temp, nrslarvae2_temp, nrslarvae3_temp) 
saveRDS(nrslarvae_COG, here('data', 'nrslarvae_COG'))

rm(df_nrslarvae1_cesm126, df_nrslarvae1_cesm585, df_nrslarvae1_gfdl126, df_nrslarvae1_gfdl585,
   df_nrslarvae1_miroc126, df_nrslarvae1_miroc585, df_nrslarvae2_cesm126, df_nrslarvae2_cesm585, 
   df_nrslarvae2_gfdl126, df_nrslarvae2_gfdl585, df_nrslarvae2_miroc126, df_nrslarvae2_miroc585,
   df_nrslarvae3_cesm126, df_nrslarvae3_cesm585, df_nrslarvae3_gfdl126, df_nrslarvae3_gfdl585,
   df_nrslarvae3_miroc126, df_nrslarvae3_miroc585, nrslarvae1_sdm, nrslarvae1_temp, nrslarvae2_sdm,
   nrslarvae2_temp, nrslarvae3_sdm, nrslarvae3_temp)

## Pacific Cod ----
### Larvae ----
# Load data
pcodlarvae_hindcast <- readRDS(here('data', 'pcod_larvae_hindcast'))
df_pcodlarvae1_cesm126 <- readRDS(here('data', 'df_pcodlarvae1_cesm126.rds'))
df_pcodlarvae2_cesm126 <- readRDS(here('data', 'df_pcodlarvae2_cesm126.rds'))
df_pcodlarvae3_cesm126 <- readRDS(here('data', 'df_pcodlarvae3_cesm126.rds'))
df_pcodlarvae1_cesm585 <- readRDS(here('data', 'df_pcodlarvae1_cesm585.rds'))
df_pcodlarvae2_cesm585 <- readRDS(here('data', 'df_pcodlarvae2_cesm585.rds'))
df_pcodlarvae3_cesm585 <- readRDS(here('data', 'df_pcodlarvae3_cesm585.rds'))
df_pcodlarvae1_gfdl126 <- readRDS(here('data', 'df_pcodlarvae1_gfdl126.rds'))
df_pcodlarvae2_gfdl126 <- readRDS(here('data', 'df_pcodlarvae2_gfdl126.rds'))
df_pcodlarvae3_gfdl126 <- readRDS(here('data', 'df_pcodlarvae3_gfdl126.rds'))
df_pcodlarvae1_gfdl585 <- readRDS(here('data', 'df_pcodlarvae1_gfdl585.rds'))
df_pcodlarvae2_gfdl585 <- readRDS(here('data', 'df_pcodlarvae2_gfdl585.rds'))
df_pcodlarvae3_gfdl585 <- readRDS(here('data', 'df_pcodlarvae3_gfdl585.rds'))
df_pcodlarvae1_miroc126 <- readRDS(here('data', 'df_pcodlarvae1_miroc126.rds'))
df_pcodlarvae2_miroc126 <- readRDS(here('data', 'df_pcodlarvae2_miroc126.rds'))
df_pcodlarvae3_miroc126 <- readRDS(here('data', 'df_pcodlarvae3_miroc126.rds'))
df_pcodlarvae1_miroc585 <- readRDS(here('data', 'df_pcodlarvae1_miroc585.rds'))
df_pcodlarvae2_miroc585 <- readRDS(here('data', 'df_pcodlarvae2_miroc585.rds'))
df_pcodlarvae3_miroc585 <- readRDS(here('data', 'df_pcodlarvae3_miroc585.rds'))

pcodlarvae1_sdm <- list(df_pcodlarvae1_cesm126[[2]], df_pcodlarvae1_cesm585[[2]],
                       df_pcodlarvae1_gfdl126[[2]], df_pcodlarvae1_gfdl585[[2]],
                       df_pcodlarvae1_miroc126[[2]], df_pcodlarvae1_miroc585[[2]])
pcodlarvae1_temp <- list(df_pcodlarvae1_cesm126[[3]], df_pcodlarvae1_cesm585[[3]],
                        df_pcodlarvae1_gfdl126[[3]], df_pcodlarvae1_gfdl585[[3]],
                        df_pcodlarvae1_miroc126[[3]], df_pcodlarvae1_miroc585[[3]])

pcodlarvae2_sdm <- list(df_pcodlarvae2_cesm126[[2]], df_pcodlarvae2_cesm585[[2]],
                       df_pcodlarvae2_gfdl126[[2]], df_pcodlarvae2_gfdl585[[2]],
                       df_pcodlarvae2_miroc126[[2]], df_pcodlarvae2_miroc585[[2]])
pcodlarvae2_temp <- list(df_pcodlarvae2_cesm126[[3]], df_pcodlarvae2_cesm585[[3]],
                        df_pcodlarvae2_gfdl126[[3]], df_pcodlarvae2_gfdl585[[3]],
                        df_pcodlarvae2_miroc126[[3]], df_pcodlarvae2_miroc585[[3]])

pcodlarvae3_sdm <- list(df_pcodlarvae3_cesm126[[2]], df_pcodlarvae3_cesm585[[2]],
                       df_pcodlarvae3_gfdl126[[2]], df_pcodlarvae3_gfdl585[[2]],
                       df_pcodlarvae3_miroc126[[2]], df_pcodlarvae3_miroc585[[2]])
pcodlarvae3_temp <- list(df_pcodlarvae3_cesm126[[3]], df_pcodlarvae3_cesm585[[3]],
                        df_pcodlarvae3_gfdl126[[3]], df_pcodlarvae3_gfdl585[[3]],
                        df_pcodlarvae3_miroc126[[3]], df_pcodlarvae3_miroc585[[3]])

pcodlarvae_COG <- COG_calc(pcodlarvae_hindcast, pcodlarvae1_sdm, pcodlarvae2_sdm, pcodlarvae3_sdm,
                          pcodlarvae1_temp, pcodlarvae2_temp, pcodlarvae3_temp) 
saveRDS(pcodlarvae_COG, here('data', 'pcodlarvae_COG'))

rm(df_pcodlarvae1_cesm126, df_pcodlarvae1_cesm585, df_pcodlarvae1_gfdl126, df_pcodlarvae1_gfdl585,
   df_pcodlarvae1_miroc126, df_pcodlarvae1_miroc585, df_pcodlarvae2_cesm126, df_pcodlarvae2_cesm585, 
   df_pcodlarvae2_gfdl126, df_pcodlarvae2_gfdl585, df_pcodlarvae2_miroc126, df_pcodlarvae2_miroc585,
   df_pcodlarvae3_cesm126, df_pcodlarvae3_cesm585, df_pcodlarvae3_gfdl126, df_pcodlarvae3_gfdl585,
   df_pcodlarvae3_miroc126, df_pcodlarvae3_miroc585, pcodlarvae1_sdm, pcodlarvae1_temp, pcodlarvae2_sdm,
   pcodlarvae2_temp, pcodlarvae3_sdm, pcodlarvae3_temp)

## Visualization and Calc ----
# Pollock
plot_COG(pkegg_COG)
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_egg_COG.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()
pkegg_COG[[10]]
mean(pkegg_COG[[9]]$distance, na.rm = TRUE)

plot_COG(pklarvae_COG)
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_COG.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

pklarvae_COG[[10]]
mean(pklarvae_COG[[9]]$distance, na.rm = TRUE)

lifestage_dist(pkegg_COG, pklarvae_COG)
temp_dist(pkegg_COG)
temp_dist(pklarvae_COG)

# Flathead Sole
plot_COG(fhsegg_COG)
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_egg_COG.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()
fhsegg_COG[[10]]
mean(fhsegg_COG[[9]]$distance, na.rm = TRUE)

plot_COG(fhslarvae_COG)
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_COG.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()

fhslarvae_COG[[10]]
mean(fhslarvae_COG[[9]]$distance, na.rm = TRUE)

lifestage_dist(fhsegg_COG, fhslarvae_COG)
temp_dist(fhsegg_COG)
temp_dist(fhslarvae_COG)

# Alaska Plaice
plot_COG(akpegg_COG)
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_egg_COG.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()
akpegg_COG[[10]]
mean(akpegg_COG[[9]]$distance, na.rm = TRUE)

plot_COG(akplarvae_COG)
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_COG.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()
akplarvae_COG[[10]]
mean(akplarvae_COG[[9]]$distance, na.rm = TRUE)

lifestage_dist(akpegg_COG, akplarvae_COG)
temp_dist(akpegg_COG)
temp_dist(akplarvae_COG)

# Yellowfin Sole
plot_COG(yfslarvae_COG)
dev.copy(jpeg,
         here('results/yellowfin_forecast',
              'yellowfin_larvae_COG.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()
yfslarvae_COG[[10]]
mean(yfslarvae_COG[[9]]$distance, na.rm = TRUE)

temp_dist(yfslarvae_COG)

# Northern Rock Sole
plot_COG(nrslarvae_COG)
dev.copy(jpeg,
         here('results/rocksole_forecast',
              'rocksole_larvae_COG.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()
nrslarvae_COG[[10]]
mean(nrslarvae_COG[[9]]$distance, na.rm = TRUE)

temp_dist(nrslarvae_COG)

# Pacific Cod
plot_COG(pcodlarvae_COG)
dev.copy(jpeg,
         here('results/cod_forecast',
              'cod_larvae_COG.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()
pcodlarvae_COG[[10]]
mean(pcodlarvae_COG[[9]]$distance, na.rm = TRUE)

temp_dist(pcodlarvae_COG)

# Velocity
# Pinsky et al. (2013) used linear regression to calculate rate of change
# Calculate COG per year, then regress latitude against year with slope being the velocity
pkegg_velocity <- velocity_calc(pkegg_COG)
pklarvae_velocity <- velocity_calc(pklarvae_COG)
fhsegg_velocity <- velocity_calc(fhsegg_COG)
fhslarvae_velocity <- velocity_calc(fhslarvae_COG)
akpegg_velocity <- velocity_calc(akpegg_COG)
akplarvae_velocity <- velocity_calc(akplarvae_COG)
yfslarvae_velocity <- velocity_calc(yfslarvae_COG)
nrslarvae_velocity <- velocity_calc(nrslarvae_COG)
pcodlarvae_velocity <- velocity_calc(pcodlarvae_COG)

hindcast <- data.frame(abundance = c(pkegg_velocity[[7]],
                                     pklarvae_velocity[[7]],
                                     fhsegg_velocity[[7]],
                                     fhslarvae_velocity[[7]],
                                     akpegg_velocity[[7]],
                                     akplarvae_velocity[[7]],
                                     yfslarvae_velocity[[7]],
                                     nrslarvae_velocity[[7]],
                                     pcodlarvae_velocity[[7]]),
                       temperature = c(pkegg_velocity[[8]],
                                       pklarvae_velocity[[8]],
                                       fhsegg_velocity[[8]],
                                       fhslarvae_velocity[[8]],
                                       akpegg_velocity[[8]],
                                       akplarvae_velocity[[8]],
                                       yfslarvae_velocity[[8]],
                                       nrslarvae_velocity[[8]],
                                       pcodlarvae_velocity[[8]]))
period1 <- data.frame(abundance = c(pkegg_velocity[[1]],
                                    pklarvae_velocity[[1]],
                                    fhsegg_velocity[[1]],
                                    fhslarvae_velocity[[1]],
                                    akpegg_velocity[[1]],
                                    akplarvae_velocity[[1]],
                                    yfslarvae_velocity[[1]],
                                    nrslarvae_velocity[[1]],
                                    pcodlarvae_velocity[[1]]),
                      temperature = c(pkegg_velocity[[4]],
                                      pklarvae_velocity[[4]],
                                      fhsegg_velocity[[4]],
                                      fhslarvae_velocity[[4]],
                                      akpegg_velocity[[4]],
                                      akplarvae_velocity[[4]],
                                      yfslarvae_velocity[[4]],
                                      nrslarvae_velocity[[4]],
                                      pcodlarvae_velocity[[4]]))
period2 <- data.frame(abundance = c(pkegg_velocity[[2]],
                                    pklarvae_velocity[[2]],
                                    fhsegg_velocity[[2]],
                                    fhslarvae_velocity[[2]],
                                    akpegg_velocity[[2]],
                                    akplarvae_velocity[[2]],
                                    yfslarvae_velocity[[2]],
                                    nrslarvae_velocity[[2]],
                                    pcodlarvae_velocity[[2]]),
                      temperature = c(pkegg_velocity[[5]],
                                      pklarvae_velocity[[5]],
                                      fhsegg_velocity[[5]],
                                      fhslarvae_velocity[[5]],
                                      akpegg_velocity[[5]],
                                      akplarvae_velocity[[5]],
                                      yfslarvae_velocity[[5]],
                                      nrslarvae_velocity[[5]],
                                      pcodlarvae_velocity[[5]]))
period3 <- data.frame(abundance = c(pkegg_velocity[[3]],
                                    pklarvae_velocity[[3]],
                                    fhsegg_velocity[[3]],
                                    fhslarvae_velocity[[3]],
                                    akpegg_velocity[[3]],
                                    akplarvae_velocity[[3]],
                                    yfslarvae_velocity[[3]],
                                    nrslarvae_velocity[[3]],
                                    pcodlarvae_velocity[[3]]),
                      temperature = c(pkegg_velocity[[6]],
                                      pklarvae_velocity[[6]],
                                      fhsegg_velocity[[6]],
                                      fhslarvae_velocity[[6]],
                                      akpegg_velocity[[6]],
                                      akplarvae_velocity[[6]],
                                      yfslarvae_velocity[[6]],
                                      nrslarvae_velocity[[6]],
                                      pcodlarvae_velocity[[6]]))

library(ggplot2)
library(gridExtra)
library(ggpmisc)
hindcast_plot <- ggplot(data = hindcast,
                        aes(x = temperature,
                            y = abundance)) +
  geom_point() +
  stat_poly_line(method = "lm",
                 se = FALSE) +
  stat_poly_eq(method = "lm")

plot1 <- ggplot(data = period1, 
                aes(x = temperature, 
                    y = abundance)) +
  geom_point() +
  stat_poly_line(method = "lm",
                 se = FALSE) +
  stat_poly_eq(method = "lm")

plot2 <- ggplot(data = period2,
                aes(x = temperature,
                    y = abundance)) +
  geom_point() +
  stat_poly_line(method = "lm",
                 se = FALSE) +
  stat_poly_eq(method = "lm")

plot3 <- ggplot(data = period3,
                aes(x = temperature,
                    y = abundance)) +
  geom_point() +
  stat_poly_line(method = "lm",
                 se = FALSE) +
  stat_poly_eq(method = "lm")

grid.arrange(hindcast_plot, plot1, plot2, plot3, nrow = 1)

# Tracking climate velocity
# Negative means that they aren't shifting as fast as climate velocity
# Something up with the AKP egg temperatures
# Not sure relative bias needed since using scaled values
hindcast$bias <- sign(hindcast$temperature) * (hindcast$abundance - hindcast$temperature)
period1$bias <- sign(period1$temperature) * (period1$abundance - period1$temperature)
period2$bias <- sign(period2$temperature) * (period2$abundance - period2$temperature)
period3$bias <- sign(period3$temperature) * (period3$abundance - period3$temperature)
# hindcast$relative_bias <- (sign(hindcast$temperature) * (hindcast$abundance - hindcast$temperature)) / hindcast$temperature
# period1$relative_bias <- (sign(period1$temperature) * (period1$abundance - period1$temperature)) / period1$temperature
# period2$relative_bias <- (sign(period2$temperature) * (period2$abundance - period2$temperature)) / period2$temperature
# period3$relative_bias <- (sign(period3$temperature) * (period3$abundance - period3$temperature)) / period3$temperature

t.test(hindcast$bias)
t.test(period1$bias)
t.test(period2$bias)
t.test(period3$bias)
# t.test(hindcast$relative_bias)
# t.test(period1$relative_bias)
# t.test(period2$relative_bias)
# t.test(period3$relative_bias)

# Look at relationship between thermal niche and abundance
summary(lm(abundance ~ temperature, data = hindcast))
summary(lm(abundance ~ temperature, data = period1))
summary(lm(abundance ~ temperature, data = period2))
summary(lm(abundance ~ temperature, data = period3))

# Sandbox
data_hindcast <- as.data.frame(cbind(lapply(akpegg_COG[[1]][[2]], "[[", 1), deparse.level = 1))
data_hindcast$latitude <- as.numeric(data_hindcast$V1)
data_hindcast <- tibble::rowid_to_column(data_hindcast, "id")
