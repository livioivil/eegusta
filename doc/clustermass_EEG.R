## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  echo = TRUE,
  warnings = FALSE,
  out.width = "500px",
  dpi=150
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----1, warning=FALSE,message=FALSE-------------------------------------------
rm(list = ls())
#remotes::install_github("angeella/ARIeeg")
#remotes::install_github("jaromilfrossard/permuco4brain")
#remotes::install_github("jaromilfrossard/permuco")
library(ARIeeg) #to compute ARI for spatial-temporal clusters
library(permuco4brain) #to compute the spatial-temporal clusters
library(plotly) # to create interactive plots (optional)
library(tidyverse) #collection of R packages (e.g., ggplot2, dplyr)
library(permuco) #to compute the temporal clusters
library(hommel) #to compute ARI for temporal clusters

## ----2------------------------------------------------------------------------
load(system.file("extdata", "data_eeg_emotion.RData", package = "ARIeeg"))

## -----------------------------------------------------------------------------
str(dati)

## ----4------------------------------------------------------------------------
m <- nrow(dati$chan_info)
chans <- dati$chan_info$electrode[c(1:(m-5))]
dati <- 
  dati %>%
  eegUtils::select_elecs(chans)%>% 
  eegUtils::select_epochs(epoch_events = c(3,5))

## ----echo = FALSE-------------------------------------------------------------
knitr::kable(head(data.frame(SUBJECT = rep(seq(20), 500*2),
                        STIMULI = rep(c("3","5"), each = 500*20),
                        TIME = rep(seq(500), 20*2))))

## -----------------------------------------------------------------------------
chans <- c("Fp1", "Fp2", "F3", "F4")

dati_sel <- dati %>%
    eegUtils::select_elecs(chans)

dati_sel$signals$subj <- rep(seq(20), 500*2)
dati_sel$signals$condition <- rep(c("3","5"), each = 500*20)
dati_sel$signals$time <- rep(seq(500), 20*2)

A<-dati_sel$signals %>% 
  pivot_longer(cols = all_of(chans)) %>%
ggplot(aes(x = time, y = value)) +
  geom_line(aes(group = condition))  +
  stat_summary(
    fun = "mean", geom = "line", alpha = 1, linewidth = 1,
    aes(color = condition),show.legend = TRUE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = .17, linetype = "dotted") +
  theme(legend.position = "bottom")+
  scale_color_manual(labels = c("Disgust", "Object"), 
                     values = c("#80bfff", "#ff8080"))+
  facet_wrap(~name)
A

## ----warning=FALSE------------------------------------------------------------
plotly::ggplotly(A)    

## -----------------------------------------------------------------------------
Fp1 <- dati %>% select(Fp1)

## -----------------------------------------------------------------------------
Fp1<-dati %>%
  eegUtils::select_elecs("Fp1")%>% 
  eegUtils::select_epochs(epoch_events = c(3,5))

signal_Fp1 <- matrix(Fp1$signals$Fp1, 
              nrow = 40,
              byrow = TRUE) #remember the structure of the signals tibble
dim(signal_Fp1)

## -----------------------------------------------------------------------------
design <- 
  data.frame(participant_id =rep(seq(20), 2),
             condition =rep(c(3,5), 20)) #as in Fp1
dim(design)

## -----------------------------------------------------------------------------
f <- signal_Fp1 ~ condition + Error(participant_id/(condition))

## -----------------------------------------------------------------------------
lm_Fp1 <- clusterlm(f,data = design)
summary(lm_Fp1)

## -----------------------------------------------------------------------------
sum(lm_Fp1$multiple_comparison$condition$uncorrected$main[c(167:212), "statistic"])

## -----------------------------------------------------------------------------
plot(lm_Fp1)

## -----------------------------------------------------------------------------
praw <- lm_Fp1$multiple_comparison$condition$uncorrected$main[,2] #extract raw pvalues


output_lm <- summary(lm_Fp1)$condition

TDP <- sapply(seq(nrow(output_lm)), function(x){ 
  
  ix <- c(output_lm$start[x]:output_lm$end[x]) #set of indices hyp of interest
  
  round(discoveries(hommel(praw), ix = ix)/length(ix),3) #apply parametric ARI
  }
)

data.frame(clustermass = output_lm$`cluster mass`,
          pmass = output_lm$`P(>mass)`,
          TDP = TDP)

## -----------------------------------------------------------------------------
signal <- array(NA, dim = c(40,500, 27))

chn <- eegUtils::channel_names(dati)[1:27]

for(i in seq(27)){
  chn_utils<-dati %>%
  eegUtils::select_elecs(chn[i])%>% 
  eegUtils::select_epochs(epoch_events = c(3,5))

signal[,,i] <- matrix(as.data.frame(chn_utils$signals)[,1], 
              nrow = 40,
              byrow = TRUE)
}
dimnames(signal) <- list(NULL, NULL, chn)


## ----8------------------------------------------------------------------------
design <- 
  data.frame(participant_id =rep(seq(20), 2),
             condition =rep(c(3,5), 20)) #as before
dim(design)

## ----fig.align="center"-------------------------------------------------------
graph <- position_to_graph(as.data.frame(dati$chan_info)[c(1:27),], 
                           name = "electrode", 
                           delta = 53,
                           x = "cart_x", 
                           y = "cart_y", 
                           z = "cart_z")
graph

## -----------------------------------------------------------------------------
plot(graph)

## ----10-----------------------------------------------------------------------
f <- signal ~ condition + Error(participant_id/(condition))

## -----------------------------------------------------------------------------
model <- permuco4brain::brainperm(formula = f,
                                  data = design,
                                  graph = graph,
                                  np = 1000,
                                  multcomp = "clustermass",
                                  return_distribution = TRUE)

## -----------------------------------------------------------------------------
print(model)

## -----------------------------------------------------------------------------
cls17 <- names(model$multiple_comparison$condition$clustermass$cluster$membership[which(as.vector(model$multiple_comparison$condition$clustermass$cluster$membership)==17)])
head(cls17,n = 10)

## -----------------------------------------------------------------------------
plot(model, samples = 300)

## -----------------------------------------------------------------------------
cls17[grepl("300", cls17)]

## -----------------------------------------------------------------------------
image(model)

## -----------------------------------------------------------------------------
ARIeeg::ARIeeg(model = model, alpha = 0.3)

## -----------------------------------------------------------------------------
alpha <- 0.3
Test <- model$multiple_comparison$condition$uncorrected$distribution
dim(Test) <- c(dim(Test)[1], dim(Test)[2]*dim(Test)[3])

pv <- matrixStats::colRanks(-abs(Test)) / nrow(Test)

lambda <- pARI::lambdaOpt(pvalues = pv, 
                          family = "simes", 
                          alpha = 0.3) #compute lambda

cvOpt = pARI::criticalVector(pvalues = pv, 
                             family = "simes", 
                             lambda = lambda, 
                             alpha = 0.3) #compute l_i

plot(sort(pv[,1]), type = "l")
for(i in seq(ncol(pv))){
  
  lines(sort(pv[,i])) #plot null distribution
  
}
lines(sort(pv[,1]), col = "red") #plot observed sorted p-values
lines(cvOpt, col = "blue") #plot l_i pARI
lines((c(1:nrow(pv))*alpha)/nrow(pv), col = "green") #plot l_i ARI


## -----------------------------------------------------------------------------
ARIeeg::ARIpermEEG(model = model, family = "simes", alpha = 0.3)

