# Check how different the values are for yes vs no of a categorical cov
# Monolix gives the beta value for that parameter, so you have to obtain  
# the value for the other group.  

#Load packages and set working directory
library(readr)
library(ggplot2)
library(dplyr)
library(readxl)
library(gridExtra)

setwd("/Users/ejm6007/Desktop/NLME/B2AcReLogit50")


#Load files 
#Document with best mdoels
B2AcRe_SelectedModels <- read_xlsx("Summary_B2AcReLogit50_ModelSelection.xlsx", 
                                   sheet = "Model Selection 1")
j=1
#Prm values
B2AcReAll<-read_csv("B2AcReLogit50AllTotal.csv")

#best models for each selecttion type
selectedrows=c(13, 5, 9, 15)
parameters<-c("delta", "log10pT0", "c", "beta", "tstart")

#keep best models
rows=selectedrows[j]
B2AcRe_SelectedModels<-B2AcRe_SelectedModels[1:rows, ]

#From the matrix of prm values, create one that contains the prm values for  
#the best models
BestModelsPrmValues<-NULL
for (i in 1:rows){
  current<-B2AcReAll[which(B2AcReAll$Run==B2AcRe_SelectedModels$Run[i]), ]
  BestModelsPrmValues<-rbind(BestModelsPrmValues, current)
}

### Check population level values of log10pT0 when there is no strong neut (pop)
### vs when there there is (beta_log10pT0_Dyn_Rebound)
### log10pT0 follows a normal dist, so log10pT0_Dyn_Rebound=log10pT0_pop+beta_log10pT0_Dyn_Rebound

beta_log10pT0_Dyn<-BestModelsPrmValues[which((grepl("_log10pT0_Dynamics_", BestModelsPrmValues$param))), c(1, 2, 10, 11)]
log10pT0_pop<-BestModelsPrmValues[which((grepl("log10pT0_pop", BestModelsPrmValues$param))), c(1, 2, 10, 11)]

#Keep log10pT0_pop for the cases where f80 is a cat cov
log10pT0_pop <- log10pT0_pop[log10pT0_pop$Run %in% beta_log10pT0_Dyn$Run, ]


#Order them in increasing order of Run to make sure that they match
beta_log10pT0_Dyn <- beta_log10pT0_Dyn[order(beta_log10pT0_Dyn$Run), ]
log10pT0_pop <- log10pT0_pop[order(log10pT0_pop$Run), ]

log10pT0_Dyn_Reb= log10pT0_pop$estimate+beta_log10pT0_Dyn$estimate

#Check the difference between the log10pT0_pop and log10pT0_f80_yes
Diff_log10pT0=round(log10pT0_Dyn_Reb-log10pT0_pop$estimate, 2)


### Find the population level value of log10beta for the monkeys in rebound 
### log10beta follows a normal dist 
beta_log10beta_Dyn<-BestModelsPrmValues[which((grepl("_log10beta_Dynamics_", BestModelsPrmValues$param))), c(1, 2, 10, 11)]
log10beta_pop<-BestModelsPrmValues[which((grepl("log10beta_pop", BestModelsPrmValues$param))), c(1, 2, 10, 11)]

#Keep log10beta_pop for the cases where f80 is a cat cov
log10beta_pop <- log10beta_pop[log10beta_pop$Run %in% beta_log10beta_Dyn$Run, ]


#Order them in increasing order of Run to make sure that they match
beta_log10beta_Dyn <- beta_log10beta_Dyn[order(beta_log10beta_Dyn$Run), ]
log10beta_pop <- log10beta_pop[order(log10beta_pop$Run), ]

log10beta_Dyn_Reb= log10beta_pop$estimate+beta_log10beta_Dyn$estimate

#Check the difference between the log10beta_pop and log10pT0_f80_yes
Diff_log10beta=round(log10beta_Dyn_Reb-log10beta_pop$estimate, 2)

n=length(log10beta_Dyn_Reb)

# Combine the data into a data frame
df <- data.frame(
  Variable = rep(c("Diff beta", "Diff pT0"), times = c(length(Diff_log10beta), length(Diff_log10pT0))),
                 Value = c(abs(Diff_log10beta), abs(Diff_log10pT0))
  )
  
ggplot(df, aes(x = Variable, y = Value, fill = Variable)) +
  geom_violin(position = "dodge", trim = FALSE) +
  ylab("log10 Values") +
  ggtitle("") +
  theme_minimal()
