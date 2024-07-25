# Projecting marine species distributions during early life stages under future climate scenarios

This is part of an NPRB funded project focused on forecasting Bering Sea groundfish distributions during their early life stages (egg and larvae). Portions of this work are found in Chapter 3 of my PhD dissertation and are published in _Fish and Fisheries_, found [here](https://onlinelibrary.wiley.com/doi/pdf/10.1111/faf.12835).

### Data
- The **EcoFOCI Groundfish Surveys** are used for egg and larval data for six groundfish species. These data are collected using paired bongo nets at predetermined stations throughout the eastern Bering Sea. Data are available upon request from EcoFOCI. Data cleaning methods can be found [here](code/EcoFOCI_cleaning.R/).
- The **Bering10K ROMS-NPZ** model output. The hindcast, historical, and forecast outputs are used for these analyses. Forecasts are simulated using the CMIP6 and are [publicly available](https://beringnpz.github.io/roms-bering-sea/B10K-dataset-docs/). Three Earth System Models (ESM) with two emission scenarios (SSP) were used for model runs and selected for these forecasts to provide sufficient range of future conditions.

### Methods
#### Bering10K Bias Correction
This process is essentially where we calibrated the Bering10K model to deal with discrepancies. The equation below shows how this was done using both the hindcasted and projected outputs. A set of reference years (1980 - 2014) were the time frame from which we calculated the means and standard deviations used in this equation. The mean projected output during the reference years was subtracted from the raw forecasted output. Then, this value was multiplied by the standard deviation of hindcasted reference years divided by the standard deviation of the projected reference years. The [calculations](code/bias_correction.R/) require significant computing power (32 GB RAM or more). 

$`T_{future',y} = \overline{T}_{hindcast,ref} + \frac{\sigma_{hindcast,ref}}{\sigma_{future,ref}} (T_{future,y} - \overline{T}_{future,ref})
`$

#### Varying Coefficient Generalized Additive Models (VGAMs)
Other [choices](code/model_exploratory.R/) for modeling were evaluated, and we ultimately settled on using VGAMs with a Tweedie distribution. VGAMs allow a coefficient to vary with change in another variable. In this application, the coefficient is either location or day of year which varies with mean temperature over the Bering Sea continental shelf. Below is an example equation, with the last term being the varying coefficient term:

$`
log(CPUE + 1) = re(year) + s(lat, lon) + s(day\ of\ year) + s(depth) + s(SST) + s(lat, lon, by = temp)
`$

- To prepare models for projecting future distributions, we parameterized using the Bering10K hindcasts [here](code/hindcast_wROMS.Rmd/).
- Once parameterized, [projections](code/projections.R/) were made for each time period, species, and life stage. Predicted abundance was averaged for each time period and for each ESM in order to get plots for each SSP.
- The mean temperature value used in the variable coefficient term was calculated from the Bering10K, illustrated [here](code/ROMS_temp_index.R/).
- Center of gravity was calculated to determine broadly how egg and larval distributions will change in the future. This provides a direction and distance of movement. These calculations were done [here](code/center_of_gravity.R).
