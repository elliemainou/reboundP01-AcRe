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
#Prm values
B2AcReAll<-read_csv("B2AcReLogit50AllTotal.csv")

# Check for entries with the first character as 0
Runs<-B2AcReAll$Run
result <- sapply(Runs, function(x) substr(x, 1, 1) == "0")

x<-which(result==FALSE)

R2AcReAllResults<-B2AcReAll[x, ]

ID<-unique(R2AcReAllResults$Run)

length(ID)

AIC<-c()
for ( i in 1:length(ID)){
  x<-which(R2AcReAllResults$Run==ID[i])
  current<-R2AcReAllResults[x, ]
  AIC[i]=current$AIC[1]
}

DeltaAIC<-AIC-min(AIC)

ggplot(data = data.frame(DeltaAIC), aes(x = DeltaAIC)) +
  geom_histogram(binwidth = 6, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "",
       x = "ΔAIC",
       y = "Frequency") +
  theme_minimal()+
  theme(axis.text = element_text(size = 18),   # Adjust the size of axis text
        axis.title = element_text(size = 16),  # Adjust the size of axis titles
        plot.title = element_text(size = 18))  # Adjust the size of plot title

