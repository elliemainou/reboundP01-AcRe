# Create pie-charts with the characteristics of the best models
# For parameters whose population value INCREASES keep the ORIGINAL COLOR. 
# For parameters whose population value DECREASES use a WASHED-OUT COLOR. 

#This is for one type for model selection. Need to repeat for all 5. 

#Load packages and set working directory
library(readr)
library(ggplot2)
library(dplyr)
library(readxl)
library(gridExtra)
library(tidyr)

setwd("/Users/ejm6007/Desktop/NLME/B2AcReLogit50")


#Load files 
#Document with best mdoels
B2AcRe_SelectedModels <- read_xlsx("Summary_B2AcReLogit50_ModelSelection.xlsx", 
                                   sheet = "Model Selection 4") ###change
j=4 ###change
#Prm values
B2AcReAll<-read_csv("B2AcReLogit50All.csv")

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

#Keep prms that are the beta prms for the covariates
betaRebound<-BestModelsPrmValues[which((grepl("Dynamics_Rebound", 
                                              BestModelsPrmValues$param))), c(1, 2, 10, 11)]

# Order the matrix based on the (param) in ascending order
betaRebound <- betaRebound[order(betaRebound$param), ]

#Keep track of each parameter
beta_delta_Rebound<-betaRebound[which((grepl("delta", betaRebound$param))), ]
beta_log10pT0_Rebound<-betaRebound[which((grepl("log10pT0", betaRebound$param))), ]
beta_c_Rebound<-betaRebound[which((grepl("_c_", betaRebound$param))), ]
beta_log10beta_Rebound<-betaRebound[which((grepl("log10beta", betaRebound$param))), ]
beta_tstart_Rebound<-betaRebound[which((grepl("tstart", betaRebound$param))), ]

#Create count matrices
name<-c("ModelSelection", "Covariate", "Prm", "Increase", "Decrease")
SummaryRebound<-as.data.frame(matrix(data=NA, nrow=length(parameters), ncol=length(name)))
colnames(SummaryRebound)<-name
SummaryRebound$ModelSelection=j
SummaryRebound$Covariate="Dynamics Rebound"
SummaryRebound$Prm=parameters

SummaryRebound$Increase[1]<-length(which(beta_delta_Rebound$estimate>0))
SummaryRebound$Increase[2]<-length(which(beta_log10pT0_Rebound$estimate>0))
SummaryRebound$Increase[3]<-length(which(beta_c_Rebound$estimate>0))
SummaryRebound$Increase[4]<-length(which(beta_log10beta_Rebound$estimate>0))
SummaryRebound$Increase[5]<-length(which(beta_tstart_Rebound$estimate>0))

SummaryRebound$Decrease[1]<-length(which(beta_delta_Rebound$estimate<0))
SummaryRebound$Decrease[2]<-length(which(beta_log10pT0_Rebound$estimate<0))
SummaryRebound$Decrease[3]<-length(which(beta_c_Rebound$estimate<0))
SummaryRebound$Decrease[4]<-length(which(beta_log10beta_Rebound$estimate<0))
SummaryRebound$Decrease[5]<-length(which(beta_tstart_Rebound$estimate<0))


#Pie charts for Dynamics

par(mfrow = c(3, 2))

# Create a list to store ggplot objects
plots_list <- list()

# Loop through each row and create a pie chart without labels and numbers
for (i in 1:nrow(SummaryRebound)) {
  count_data <- data.frame(
    Change = c("Increase", "Decrease"),
    Count = c(
      SummaryRebound$Increase[i],
      SummaryRebound$Decrease[i]
    )
  )
  
  # Create and store the ggplot object in the list
  plots_list[[i]] <- ggplot(count_data, aes(x = "", y = Count, fill = Change)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_brewer(palette = "Set1") +
    theme(
      legend.position = "none",  # Remove legend
      legend.title = element_blank(),  # Remove legend title
      legend.text = element_text(size = 14)  # Adjust the size of legend text
    ) +
    ggtitle(paste("Dynamics, Selection ", j, ", ", SummaryRebound$Prm[i], sep = ""))
}

# Combine and arrange the plots side by side
grid.arrange(grobs = plots_list, ncol = 2)

# Reset the layout
par(mfrow = c(1, 1))

SummaryRebound$Increase=SummaryRebound$Increase/selectedrows[j]
SummaryRebound$Decrease=SummaryRebound$Decrease/selectedrows[j]


### Create barplot containing the proportion of best models a parameter appears
### separated by how muych it changes 
summary_long <- gather(SummaryRebound, key = "Change", value = "Value", Increase, Decrease)

# Creating a stacked bar plot
ggplot(summary_long, aes(x = Prm, y = Value, fill = Change)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = selectedrows[j], linetype = "dashed", color = "black") +
  labs(title = "",
       x = "Parameter",
       y = "Proportion of selected models") +
  scale_fill_manual(values = c("Increase" = "deepskyblue4", "Decrease" = "brown1")) +
  theme_minimal() +
  ylim(0, 1) + 
  theme(axis.text = element_text(size = 18, family = "Times New Roman"),  # Adjust the size and font for axis text
        axis.title = element_text(size = 20, family = "Times New Roman"),  # Set the font for axes titles
        plot.title = element_text(size = 24, face = "bold", family = "Times New Roman"), 
        legend.text = element_text(size = 16),  # Adjust legend text size
        legend.title = element_text(size = 18), 
        text = element_text(family = "Times New Roman"))+  # Set the font for the plot title
  scale_x_discrete(labels = c("β", "c", "δ", expression(pT[~"0"]), expression(t[~"start"])))




