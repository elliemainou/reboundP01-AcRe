#Performed a paired test to see if the values between acute and rebound are 
#different categorized by presence of absence of nAbs
# Histogram of c for B2AcReLogit50

library(deSolve)
library(ggplot2)
library(readr)
library(stringr) 
library(patchwork)
library(tidyverse)
setwd("/Users/ejm6007/Desktop/NLME/B2AcReLogit50")

PVL <- read_csv("~/Desktop/NLME/PVL_B2C5C4_ATI_f_LR.csv")
B2AcReLogit50IndivPrmAll <- read_csv("B2AcReLogit50IndivPrmAll.csv")
Summary_B2AcReLogit50 <- read_xlsx("Summary_B2AcReLogit50_ModelSelection.xlsx", 
                                   sheet = "Model Selection 1") ###change

j=1 ###change

selectedrows=c(13, 5, 9, 15)

#Find indiv prms of best models 

ID<-Summary_B2AcReLogit50$Run[1:selectedrows[j]]


BestModels<-NULL
for(i in 1:length(ID)){
  current<-B2AcReLogit50IndivPrmAll[which(B2AcReLogit50IndivPrmAll$x==ID[i]), ]
  BestModels<-rbind(BestModels,current)
}

Run<-unique(BestModels$x)

for (j in 1:nrow(BestModels)){
  strlength<-str_length(BestModels$id[j])
  Dynamics_current<-str_sub(BestModels$id[j], 7,strlength)
  BestModels$Dynamics[j]<-Dynamics_current
}

#Add which monkeys have f80
for (j in 1:nrow(BestModels)){
  id_current<-str_sub(BestModels$id[j], 1,5)
  f80=PVL$f80[which(PVL$animal_id==id_current)[1]]
  BestModels$f80[j]<-f80
}

for (i in 1:length(Run)){
  current <- BestModels[which(BestModels$x == Run[i]), ]
  data_yes <- subset(current, f80 == "yes")
  
  
  wilcoxon_test_result_yes <- wilcox.test(data_yes$delta[which(data_yes$Dynamics=="Rebound")], 
                                          data_yes$delta[which(data_yes$Dynamics=="Acute")], paired = TRUE)
  
  print(wilcoxon_test_result_yes)
  
  test.result[i, 1]<-wilcoxon_test_result_yes$p.value
}

for (i in 1:length(Run)){
  current <- BestModels[which(BestModels$x == Run[i]), ]
  data_yes <- subset(current, f80 == "no")
  
  
  wilcoxon_test_result_yes <- wilcox.test(data_yes$delta[which(data_yes$Dynamics=="Rebound")], 
                                          data_yes$delta[which(data_yes$Dynamics=="Acute")], paired = TRUE)
  
  print(wilcoxon_test_result_yes)
  
  test.result[i, 2]<-wilcoxon_test_result_yes$p.value
}
