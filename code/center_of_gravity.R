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

pkegg1_list <- list(df_pkegg1_cesm126, df_pkegg1_cesm585,
                   df_pkegg1_gfdl126, df_pkegg1_gfdl585,
                   df_pkegg1_miroc126, df_pkegg1_miroc585)

pkegg2_list <- list(df_pkegg2_cesm126, df_pkegg2_cesm585,
                   df_pkegg2_gfdl126, df_pkegg2_gfdl585,
                   df_pkegg2_miroc126, df_pkegg2_miroc585)

pkegg3_list <- list(df_pkegg3_cesm126, df_pkegg3_cesm585,
                   df_pkegg3_gfdl126, df_pkegg3_gfdl585,
                   df_pkegg3_miroc126, df_pkegg3_miroc585)

pkegg_COG <- COG_time(pkegg_hindcast, pkegg1_list, pkegg2_list, pkegg3_list)
saveRDS(pkegg_COG, here('data', 'pkegg_COG'))

rm(df_pkegg1_cesm126, df_pkegg1_cesm585, df_pkegg1_gfdl126, df_pkegg1_gfdl585,
   df_pkegg1_miroc126, df_pkegg1_miroc585, df_pkegg2_cesm126, df_pkegg2_cesm585, 
   df_pkegg2_gfdl126, df_pkegg2_gfdl585, df_pkegg2_miroc126, df_pkegg2_miroc585,
   df_pkegg3_cesm126, df_pkegg3_cesm585, df_pkegg3_gfdl126, df_pkegg3_gfdl585,
   df_pkegg3_miroc126, df_pkegg3_miroc585, pkegg1_list, pkegg2_list, pkegg3_list,
   pkegg_hindcast)

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

pklarvae1_list <- list(df_pklarvae1_cesm126, df_pklarvae1_cesm585,
                    df_pklarvae1_gfdl126, df_pklarvae1_gfdl585,
                    df_pklarvae1_miroc126, df_pklarvae1_miroc585)

pklarvae2_list <- list(df_pklarvae2_cesm126, df_pklarvae2_cesm585,
                    df_pklarvae2_gfdl126, df_pklarvae2_gfdl585,
                    df_pklarvae2_miroc126, df_pklarvae2_miroc585)

pklarvae3_list <- list(df_pklarvae3_cesm126, df_pklarvae3_cesm585,
                    df_pklarvae3_gfdl126, df_pklarvae3_gfdl585,
                    df_pklarvae3_miroc126, df_pklarvae3_miroc585)

pklarvae_COG <- COG_time(pklarvae_hindcast, pklarvae1_list, pklarvae2_list, pklarvae3_list)
saveRDS(pklarvae_COG, here('data', 'pklarvae_COG'))

rm(df_pklarvae1_cesm126, df_pklarvae1_cesm585, df_pklarvae1_gfdl126, df_pklarvae1_gfdl585,
   df_pklarvae1_miroc126, df_pklarvae1_miroc585, df_pklarvae2_cesm126, df_pklarvae2_cesm585, 
   df_pklarvae2_gfdl126, df_pklarvae2_gfdl585, df_pklarvae2_miroc126, df_pklarvae2_miroc585,
   df_pklarvae3_cesm126, df_pklarvae3_cesm585, df_pklarvae3_gfdl126, df_pklarvae3_gfdl585,
   df_pklarvae3_miroc126, df_pklarvae3_miroc585, pklarvae1_list, pklarvae2_list, pklarvae3_list,
   pklarvae_hindcast)

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

fhsegg1_list <- list(df_fhsegg1_cesm126, df_fhsegg1_cesm585,
                    df_fhsegg1_gfdl126, df_fhsegg1_gfdl585,
                    df_fhsegg1_miroc126, df_fhsegg1_miroc585)

fhsegg2_list <- list(df_fhsegg2_cesm126, df_fhsegg2_cesm585,
                    df_fhsegg2_gfdl126, df_fhsegg2_gfdl585,
                    df_fhsegg2_miroc126, df_fhsegg2_miroc585)

fhsegg3_list <- list(df_fhsegg3_cesm126, df_fhsegg3_cesm585,
                    df_fhsegg3_gfdl126, df_fhsegg3_gfdl585,
                    df_fhsegg3_miroc126, df_fhsegg3_miroc585)

fhsegg_COG <- COG_time(fhsegg_hindcast, fhsegg1_list, fhsegg2_list, fhsegg3_list)
saveRDS(fhsegg_COG, here('data', 'fhsegg_COG'))

rm(df_fhsegg1_cesm126, df_fhsegg1_cesm585, df_fhsegg1_gfdl126, df_fhsegg1_gfdl585,
   df_fhsegg1_miroc126, df_fhsegg1_miroc585, df_fhsegg2_cesm126, df_fhsegg2_cesm585, 
   df_fhsegg2_gfdl126, df_fhsegg2_gfdl585, df_fhsegg2_miroc126, df_fhsegg2_miroc585,
   df_fhsegg3_cesm126, df_fhsegg3_cesm585, df_fhsegg3_gfdl126, df_fhsegg3_gfdl585,
   df_fhsegg3_miroc126, df_fhsegg3_miroc585, fhsegg1_list, fhsegg2_list, fhsegg3_list,
   fhsegg_hindcast)

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

fhslarvae1_list <- list(df_fhslarvae1_cesm126, df_fhslarvae1_cesm585,
                       df_fhslarvae1_gfdl126, df_fhslarvae1_gfdl585,
                       df_fhslarvae1_miroc126, df_fhslarvae1_miroc585)

fhslarvae2_list <- list(df_fhslarvae2_cesm126, df_fhslarvae2_cesm585,
                       df_fhslarvae2_gfdl126, df_fhslarvae2_gfdl585,
                       df_fhslarvae2_miroc126, df_fhslarvae2_miroc585)

fhslarvae3_list <- list(df_fhslarvae3_cesm126, df_fhslarvae3_cesm585,
                       df_fhslarvae3_gfdl126, df_fhslarvae3_gfdl585,
                       df_fhslarvae3_miroc126, df_fhslarvae3_miroc585)

fhslarvae_COG <- COG_time(fhslarvae_hindcast, fhslarvae1_list, fhslarvae2_list, fhslarvae3_list)
saveRDS(fhslarvae_COG, here('data', 'fhslarvae_COG'))

rm(df_fhslarvae1_cesm126, df_fhslarvae1_cesm585, df_fhslarvae1_gfdl126, df_fhslarvae1_gfdl585,
   df_fhslarvae1_miroc126, df_fhslarvae1_miroc585, df_fhslarvae2_cesm126, df_fhslarvae2_cesm585, 
   df_fhslarvae2_gfdl126, df_fhslarvae2_gfdl585, df_fhslarvae2_miroc126, df_fhslarvae2_miroc585,
   df_fhslarvae3_cesm126, df_fhslarvae3_cesm585, df_fhslarvae3_gfdl126, df_fhslarvae3_gfdl585,
   df_fhslarvae3_miroc126, df_fhslarvae3_miroc585, fhslarvae1_list, fhslarvae2_list, fhslarvae3_list,
   fhslarvae_hindcast)

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

akpegg1_list <- list(df_akpegg1_cesm126, df_akpegg1_cesm585,
                     df_akpegg1_gfdl126, df_akpegg1_gfdl585,
                     df_akpegg1_miroc126, df_akpegg1_miroc585)

akpegg2_list <- list(df_akpegg2_cesm126, df_akpegg2_cesm585,
                     df_akpegg2_gfdl126, df_akpegg2_gfdl585,
                     df_akpegg2_miroc126, df_akpegg2_miroc585)

akpegg3_list <- list(df_akpegg3_cesm126, df_akpegg3_cesm585,
                     df_akpegg3_gfdl126, df_akpegg3_gfdl585,
                     df_akpegg3_miroc126, df_akpegg3_miroc585)

akpegg_COG <- COG_time(akpegg_hindcast, akpegg1_list, akpegg2_list, akpegg3_list)
saveRDS(akpegg_COG, here('data', 'akpegg_COG'))

rm(df_akpegg1_cesm126, df_akpegg1_cesm585, df_akpegg1_gfdl126, df_akpegg1_gfdl585,
   df_akpegg1_miroc126, df_akpegg1_miroc585, df_akpegg2_cesm126, df_akpegg2_cesm585, 
   df_akpegg2_gfdl126, df_akpegg2_gfdl585, df_akpegg2_miroc126, df_akpegg2_miroc585,
   df_akpegg3_cesm126, df_akpegg3_cesm585, df_akpegg3_gfdl126, df_akpegg3_gfdl585,
   df_akpegg3_miroc126, df_akpegg3_miroc585, akpegg1_list, akpegg2_list, akpegg3_list,
   akpegg_hindcast)

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

akplarvae1_list <- list(df_akplarvae1_cesm126, df_akplarvae1_cesm585,
                        df_akplarvae1_gfdl126, df_akplarvae1_gfdl585,
                        df_akplarvae1_miroc126, df_akplarvae1_miroc585)

akplarvae2_list <- list(df_akplarvae2_cesm126, df_akplarvae2_cesm585,
                        df_akplarvae2_gfdl126, df_akplarvae2_gfdl585,
                        df_akplarvae2_miroc126, df_akplarvae2_miroc585)

akplarvae3_list <- list(df_akplarvae3_cesm126, df_akplarvae3_cesm585,
                        df_akplarvae3_gfdl126, df_akplarvae3_gfdl585,
                        df_akplarvae3_miroc126, df_akplarvae3_miroc585)

akplarvae_COG <- COG_time(akplarvae_hindcast, akplarvae1_list, akplarvae2_list, akplarvae3_list)
saveRDS(akplarvae_COG, here('data', 'akplarvae_COG'))

rm(df_akplarvae1_cesm126, df_akplarvae1_cesm585, df_akplarvae1_gfdl126, df_akplarvae1_gfdl585,
   df_akplarvae1_miroc126, df_akplarvae1_miroc585, df_akplarvae2_cesm126, df_akplarvae2_cesm585, 
   df_akplarvae2_gfdl126, df_akplarvae2_gfdl585, df_akplarvae2_miroc126, df_akplarvae2_miroc585,
   df_akplarvae3_cesm126, df_akplarvae3_cesm585, df_akplarvae3_gfdl126, df_akplarvae3_gfdl585,
   df_akplarvae3_miroc126, df_akplarvae3_miroc585, akplarvae1_list, akplarvae2_list, akplarvae3_list,
   akplarvae_hindcast)

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

yfslarvae1_list <- list(df_yfslarvae1_cesm126, df_yfslarvae1_cesm585,
                        df_yfslarvae1_gfdl126, df_yfslarvae1_gfdl585,
                        df_yfslarvae1_miroc126, df_yfslarvae1_miroc585)

yfslarvae2_list <- list(df_yfslarvae2_cesm126, df_yfslarvae2_cesm585,
                        df_yfslarvae2_gfdl126, df_yfslarvae2_gfdl585,
                        df_yfslarvae2_miroc126, df_yfslarvae2_miroc585)

yfslarvae3_list <- list(df_yfslarvae3_cesm126, df_yfslarvae3_cesm585,
                        df_yfslarvae3_gfdl126, df_yfslarvae3_gfdl585,
                        df_yfslarvae3_miroc126, df_yfslarvae3_miroc585)

yfslarvae_COG <- COG_time(yfslarvae_hindcast, yfslarvae1_list, yfslarvae2_list, yfslarvae3_list)
saveRDS(yfslarvae_COG, here('data', 'yfslarvae_COG'))

rm(df_yfslarvae1_cesm126, df_yfslarvae1_cesm585, df_yfslarvae1_gfdl126, df_yfslarvae1_gfdl585,
   df_yfslarvae1_miroc126, df_yfslarvae1_miroc585, df_yfslarvae2_cesm126, df_yfslarvae2_cesm585, 
   df_yfslarvae2_gfdl126, df_yfslarvae2_gfdl585, df_yfslarvae2_miroc126, df_yfslarvae2_miroc585,
   df_yfslarvae3_cesm126, df_yfslarvae3_cesm585, df_yfslarvae3_gfdl126, df_yfslarvae3_gfdl585,
   df_yfslarvae3_miroc126, df_yfslarvae3_miroc585, yfslarvae1_list, yfslarvae2_list, yfslarvae3_list,
   yfslarvae_hindcast)

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

nrslarvae1_list <- list(df_nrslarvae1_cesm126, df_nrslarvae1_cesm585,
                        df_nrslarvae1_gfdl126, df_nrslarvae1_gfdl585,
                        df_nrslarvae1_miroc126, df_nrslarvae1_miroc585)

nrslarvae2_list <- list(df_nrslarvae2_cesm126, df_nrslarvae2_cesm585,
                        df_nrslarvae2_gfdl126, df_nrslarvae2_gfdl585,
                        df_nrslarvae2_miroc126, df_nrslarvae2_miroc585)

nrslarvae3_list <- list(df_nrslarvae3_cesm126, df_nrslarvae3_cesm585,
                        df_nrslarvae3_gfdl126, df_nrslarvae3_gfdl585,
                        df_nrslarvae3_miroc126, df_nrslarvae3_miroc585)

nrslarvae_COG <- COG_time(nrslarvae_hindcast, nrslarvae1_list, nrslarvae2_list, nrslarvae3_list)
saveRDS(nrslarvae_COG, here('data', 'nrslarvae_COG'))

rm(df_nrslarvae1_cesm126, df_nrslarvae1_cesm585, df_nrslarvae1_gfdl126, df_nrslarvae1_gfdl585,
   df_nrslarvae1_miroc126, df_nrslarvae1_miroc585, df_nrslarvae2_cesm126, df_nrslarvae2_cesm585, 
   df_nrslarvae2_gfdl126, df_nrslarvae2_gfdl585, df_nrslarvae2_miroc126, df_nrslarvae2_miroc585,
   df_nrslarvae3_cesm126, df_nrslarvae3_cesm585, df_nrslarvae3_gfdl126, df_nrslarvae3_gfdl585,
   df_nrslarvae3_miroc126, df_nrslarvae3_miroc585, nrslarvae1_list, nrslarvae2_list, nrslarvae3_list,
   nrslarvae_hindcast)

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

pcodlarvae1_list <- list(df_pcodlarvae1_cesm126, df_pcodlarvae1_cesm585,
                        df_pcodlarvae1_gfdl126, df_pcodlarvae1_gfdl585,
                        df_pcodlarvae1_miroc126, df_pcodlarvae1_miroc585)

pcodlarvae2_list <- list(df_pcodlarvae2_cesm126, df_pcodlarvae2_cesm585,
                        df_pcodlarvae2_gfdl126, df_pcodlarvae2_gfdl585,
                        df_pcodlarvae2_miroc126, df_pcodlarvae2_miroc585)

pcodlarvae3_list <- list(df_pcodlarvae3_cesm126, df_pcodlarvae3_cesm585,
                        df_pcodlarvae3_gfdl126, df_pcodlarvae3_gfdl585,
                        df_pcodlarvae3_miroc126, df_pcodlarvae3_miroc585)

pcodlarvae_COG <- COG_time(pcodlarvae_hindcast, pcodlarvae1_list, pcodlarvae2_list, pcodlarvae3_list)
saveRDS(pcodlarvae_COG, here('data', 'pcodlarvae_COG'))

rm(df_pcodlarvae1_cesm126, df_pcodlarvae1_cesm585, df_pcodlarvae1_gfdl126, df_pcodlarvae1_gfdl585,
   df_pcodlarvae1_miroc126, df_pcodlarvae1_miroc585, df_pcodlarvae2_cesm126, df_pcodlarvae2_cesm585, 
   df_pcodlarvae2_gfdl126, df_pcodlarvae2_gfdl585, df_pcodlarvae2_miroc126, df_pcodlarvae2_miroc585,
   df_pcodlarvae3_cesm126, df_pcodlarvae3_cesm585, df_pcodlarvae3_gfdl126, df_pcodlarvae3_gfdl585,
   df_pcodlarvae3_miroc126, df_pcodlarvae3_miroc585, pcodlarvae1_list, pcodlarvae2_list, pcodlarvae3_list,
   pcodlarvae_hindcast)

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
pkegg_COG[[2]]
mean(pkegg_COG[[1]]$distance, na.rm = TRUE)

plot_COG(pklarvae_COG)
dev.copy(jpeg,
         here('results/pollock_forecast',
              'pollock_larvae_COG.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()
pklarvae_COG[[2]]
mean(pklarvae_COG[[1]]$distance, na.rm = TRUE)

lifestage_dist(pkegg_COG, pklarvae_COG)

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
fhsegg_COG[[2]]
mean(fhsegg_COG[[1]]$distance, na.rm = TRUE)

plot_COG(fhslarvae_COG)
dev.copy(jpeg,
         here('results/flathead_forecast',
              'flathead_larvae_COG.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()
fhslarvae_COG[[2]]
mean(fhslarvae_COG[[1]]$distance, na.rm = TRUE)

lifestage_dist(fhsegg_COG, fhslarvae_COG)

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
akpegg_COG[[2]]
mean(akpegg_COG[[1]]$distance, na.rm = TRUE)

plot_COG(akplarvae_COG)
dev.copy(jpeg,
         here('results/plaice_forecast',
              'plaice_larvae_COG.jpg'),
         height = 6,
         width = 6,
         res = 200,
         units = 'in')
dev.off()
akplarvae_COG[[2]]
mean(akplarvae_COG[[1]]$distance, na.rm = TRUE)

lifestage_dist(akpegg_COG, akplarvae_COG)

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
yfslarvae_COG[[2]]
mean(yfslarvae_COG[[1]]$distance, na.rm = TRUE)

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
nrslarvae_COG[[2]]
mean(nrslarvae_COG[[1]]$distance, na.rm = TRUE)

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
pcodlarvae_COG[[2]]
mean(pcodlarvae_COG[[1]]$distance, na.rm = TRUE)