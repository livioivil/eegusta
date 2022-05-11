# remotes::install_github("bnicenboim/eeguana@experimental")

rm(list=ls())
setwd("~/hdeeg")
# remotes::install_github("bnicenboim/eeguana@experimental")
library(dplyr)
library(eeguana)
library(parallel)
devtools::install("~/eegusta/")
library(eegusta)

files=dir("./data_bdf",full.names = TRUE,pattern = "bdf$")

# preproc_subj(files[14])
mclapply(1:length(files), function(i)preproc_subj(files[i],overwrite = TRUE), mc.cores = min(7,length(files)))

