## Create violin plots for the difference between acute and rebound


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
                                   sheet = "Model Selection 1") ###change
j=1 ###change
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
beta_log10pT0_Rebound<-betaRebound[which((grepl("log10pT0", betaRebound$param))), ]
beta_log10beta_Rebound<-betaRebound[which((grepl("log10beta", betaRebound$param))), ]

log10pT0_Rebound<-NULL
log10pT0_Acute<-NULL
for (i in 1:nrow(beta_log10pT0_Rebound)){
  current<-BestModelsPrmValues[which(BestModelsPrmValues$Run==beta_log10pT0_Rebound$Run[i]), ]
  
  log10pT0_pop<-current$estimate[which(current$param=="log10pT0_pop")]
  log10pT0_Acute<-rbind(log10pT0_Acute, log10pT0_pop)
  
  log10pT0<-current$estimate[which(current$param=="log10pT0_pop")]+
    current$estimate[which(current$param=="beta_log10pT0_Dynamics_Rebound")]
  log10pT0_Rebound<-rbind(log10pT0_Rebound, log10pT0)
}

log10pT0_diff<-abs(log10pT0_Rebound-log10pT0_Acute)


log10beta_Rebound<-NULL
log10beta_Acute<-NULL
for (i in 1:nrow(beta_log10beta_Rebound)){
  current<-BestModelsPrmValues[which(BestModelsPrmValues$Run==beta_log10beta_Rebound$Run[i]), ]
  
  log10beta_pop<-current$estimate[which(current$param=="log10beta_pop")]
  log10beta_Acute<-rbind(log10beta_Acute, log10beta_pop)
  
  log10beta<-current$estimate[which(current$param=="log10beta_pop")]+
    current$estimate[which(current$param=="beta_log10beta_Dynamics_Rebound")]
  log10beta_Rebound<-rbind(log10beta_Rebound, log10beta)
}

log10beta_diff<-abs(log10beta_Rebound-log10beta_Acute)



# Create a data frame
data <- data.frame(
  Parameter = rep(c("pT0", "β"), 
                  times = c(length(log10pT0_diff), length(log10beta_diff))),
  Value = c(log10pT0_diff, log10beta_diff)
)

# Create violin plot
ggplot(data, aes(x = Parameter, y = Value, fill = Parameter)) +
  geom_violin() +
  theme_minimal() +
  scale_fill_manual(values = c("brown1", "deepskyblue4")) +
  labs(title = "",
       x = "Parameter",
       y = expression("Absolute value " ~ log[10]("Rebound") - log[10]("Acute")))+ 
  theme(axis.text = element_text(size = 18, family = "Times New Roman"),  # Adjust the size and font for axis text
        axis.title = element_text(size = 20, family = "Times New Roman"),  # Set the font for axes titles
        plot.title = element_text(size = 24, face = "bold", family = "Times New Roman"), 
        legend.text = element_text(size = 16),  # Adjust legend text size
        legend.title = element_text(size = 18), 
        text = element_text(family = "Times New Roman"))+  # Set the font for the plot title
  scale_x_discrete(labels = c("pT0" = expression(pT[0]), "β" = "β"))



#Make scatterplots 
# Assuming log10pT0_Acute, log10pT0_Rebound, log10beta_Acute, and log10beta_Rebound are your vectors
result_pT0 <- t.test(log10pT0_Acute, log10pT0_Rebound)
result_beta <- t.test(log10beta_Acute, log10beta_Rebound)

# Create data frames for plotting
data_pT0 <- data.frame(
  Value = c(log10pT0_Acute, log10pT0_Rebound),
  Parameter = rep(c("Acute", "Rebound"), each = length(log10pT0_Acute)),
  Color = rep(c("brown1", "deepskyblue4"), each = length(log10pT0_Acute)),
  Variable = rep("log10pT0", length(log10pT0_Acute) * 2)
)

data_beta <- data.frame(
  Value = c(log10beta_Acute, log10beta_Rebound),
  Parameter = rep(c("Acute", "Rebound"), each = length(log10beta_Acute)),
  Color = rep(c("brown1", "deepskyblue4"), each = length(log10beta_Acute)),
  Variable = rep("log10beta", length(log10beta_Acute) * 2)
)


# Plot log10pT0 side by side with annotation
plot_pT0 <- ggplot(data_pT0, aes(x = Parameter, y = Value, color = Color)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  theme_minimal() +
  labs(title = expression(Log[10](pT[0])),
       x = "Infection Stage",
       y = expression("Mean value " (day^{-1}))) +
  theme(axis.text = element_text(size = 18, family = "Times New Roman"),  # Adjust the size and font for axis text
        axis.title = element_text(size = 20, family = "Times New Roman"),  # Set the font for axes titles
        plot.title = element_text(size = 22, face = "bold", family = "Times New Roman"), 
        legend.text = element_text(size = 16),  # Adjust legend text size
        legend.title = element_text(size = 18), 
        text = element_text(family = "Times New Roman"))+  # Set the font for the plot title
  scale_color_identity() +  # Use the specified colors without any scaling
  annotate("text", x = 1.5, y = max(data_pT0$Value) + 1, 
           label = paste("p-value =", signif(result_pT0$p.value, digits = 3)), 
           size = 5, color = "black", family = "Times New Roman") + # Add p-value annotation
  scale_x_discrete(labels = c("Acute" = "Acute infection", "Rebound" = "Viral rebound"))

# Plot log10beta side by side with annotation
plot_beta <- ggplot(data_beta, aes(x = Parameter, y = Value, color = Color)) +
  geom_point(position = position_dodge(width = 0.8), size = 3) +
  theme_minimal() +
  labs(title = expression(Log[10]("β")),
       x = "Infection Stage",
       y = "Mean value (mL/day)") +
  theme(axis.text = element_text(size = 18, family = "Times New Roman"),  # Adjust the size and font for axis text
        axis.title = element_text(size = 20, family = "Times New Roman"),  # Set the font for axes titles
        plot.title = element_text(size = 22, face = "bold", family = "Times New Roman"), 
        legend.text = element_text(size = 16),  # Adjust legend text size
        legend.title = element_text(size = 18), 
        text = element_text(family = "Times New Roman"))+  # Set the font for the plot title
  scale_color_identity() +  # Use the specified colors without any scaling
  annotate("text", x = 1.5, y = max(data_beta$Value) + 1, 
           label = paste("p-value =", signif(result_beta$p.value, digits = 3)), 
           size = 5, color = "black", family = "Times New Roman") + # Add p-value annotation
scale_x_discrete(labels = c("Acute" = "Acute infection", "Rebound" = "Viral rebound"))

# Arrange the plots side by side
grid.arrange(plot_pT0, plot_beta, ncol = 2)


