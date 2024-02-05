####################################################################################################
# PACKAGES
library(plyr)
library(dplyr)
library(readxl)
library(xml2)
library(stringr)
library(Rmisc)
library(doBy)
library(tidyverse)
library(tidyr)
library(lubridate)

library(R.matlab)

library(lme4)
library(car)
library(lmerTest)
library(pbkrtest)
library(emmeans)
library(BayesFactor)

library(colorspace)
library(ggplot2)
library(cowplot)
library(ggnewscale)
library(visreg)
library(gghalves)
library(ggpubr)
library(ggsci)
library(pals)
library(ggridges)

library(knitr)
library(pander)
library(huxtable)
library(flextable)

# set debugging options like a normal person....
options(error = function() {
    sink(stderr())
    on.exit(sink(NULL))
    traceback(3, max.lines = 1L)
    if (!interactive()) {
        q(status = 1)
    }
})
####################################################################################################
# PATHS 

host_name <- Sys.info()[4]
if (host_name == 'tux'){
    experiment_path <- '/datadisk/projects_backed_up/autocorrelation'
} else if (host_name == 'tomo-office-desktop'){
    experiment_path <-  '/DATA2/autocorrelation'
} else if (host_name == 'mac-BX23-117.local'){
    experiment_path <-  '~/projects_backed_up/acf'
}

####################################################################################################
# PARAMS 

z_snr_alpha = 0.01;

fontsize <- 14

col_eeg <- '#4000a6'


