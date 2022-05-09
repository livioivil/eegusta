# remotes::install_github("bnicenboim/eeguana@experimental")

rm(list=ls())
setwd("~/hdeeg")
setwd('C:/Users/livio/Dropbox (unipd)/dataset/hdeeg_emotion')

# remotes::install_github("bnicenboim/eeguana@experimental")
library(dplyr)
library(eeguana)

range_s=c(.05, .15)
file_name="data_seg_median_sample.Rdata"

file_out=strsplit(file_name,"\\.")[[1]]
file_out= paste0(file_out[1:(length(file_out)-1)],"_avepot",paste0(range_s[1]*100,"_",range_s[2]*100),".Rdata")


load(file_name)

D=data_seg%>%filter(between(as_time(.sample, unit = "s"), range_s[1], range_s[2])) %>%
  group_by(condition,subj)%>%summarize_all(mean, na.rm = TRUE) 
dim(D$.signal)
D$.signal$.id
D$.segments
table(D$.segments$subj,D$.segments$condition)


save(file=file_out,D)


