library(readr)
library(ggplot2)
library(dplyr)
library(readxl) 
library(RColorBrewer)

setwd("/Users/ejm6007/Desktop/NLME/B2AcReLogit50")

Summary_B2AcReLogit50Total<-read_excel("Summary_B2AcReLogit50_ModelSelection.xlsx", 
                       sheet = "Model Selection 4") 
j=4
selectedrows=c(13, 5, 9, 15)
rows=selectedrows[j] #1:rows are the best models 
Summary_B2AcReLogit50Total<-Summary_B2AcReLogit50Total[1:rows, ]

#Histogram of c values
cvalues <- Summary_B2AcReLogit50Total$cvalue
# Convert 'cvalues' to a factor with the desired order
cvalues <- factor(cvalues, levels = c(3, 6, 12))

# Create a barplot
ggplot(data.frame(cvalues), aes(x = cvalues)) +
  geom_bar(fill = "blue3") +
  labs(title = "", x = "c", y = "Count") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),    # Adjust the size of axis labels
    axis.title = element_text(size = 20)    # Adjust the size of axis titles
  )

## Create barplot when c differs between groups
yesc<-Summary_B2AcReLogit50Total[which(Summary_B2AcReLogit50Total$"c"=="Yes"), ]

cvalues <- yesc$cvalue
# Convert 'cvalues' to a factor with the desired order
cvalues <- factor(cvalues, levels = c(3, 6, 12))

# Create a barplot
ggplot(data.frame(cvalues), aes(x = cvalues)) +
  geom_bar(fill = "blue3") +
  labs(title = "", x = "c for acute infection", y = "Count") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),    # Adjust the size of axis labels
    axis.title = element_text(size = 20)    # Adjust the size of axis titles
  )


## Create barplot when c doesn't differs between groups
noc<-Summary_B2AcReLogit50Total[which(is.na(Summary_B2AcReLogit50Total$"c...13")), ]

cvalues <- noc$cvalue
# Convert 'cvalues' to a factor with the desired order
cvalues <- factor(cvalues, levels = c(3, 6, 12))

# Create a barplot
ggplot(data.frame(cvalues), aes(x = cvalues)) +
  geom_bar(fill = "blue3") +
  labs(title = "", x = "c", y = "Count") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),    # Adjust the size of axis labels
    axis.title = element_text(size = 20)    # Adjust the size of axis titles
  )



#pie chart for Dynamics

delta_Dyn<-which(Summary_B2AcReLogit50Total$delta=="Yes")
log10pT0_Dyn<-which(Summary_B2AcReLogit50Total$log10pT0=="Yes")
c_Dyn<-which(Summary_B2AcReLogit50Total$"c"=="Yes")
log10beta_Dyn<-which(Summary_B2AcReLogit50Total$log10beta=="Yes")
tstart_Dyn<-which(Summary_B2AcReLogit50Total$tstart=="Yes")

total_Dyn<-sum(delta_Dyn, log10pT0_Dyn, c_Dyn, log10beta_Dyn, tstart_Dyn)

# Create a data frame with the counts
count_data <- data.frame(
  Variable = c("delta", "pT0", "c", "beta", "tstart"),
  Count = c(
    length(delta_Dyn),
    length(log10pT0_Dyn),
    length(c_Dyn),
    length(log10beta_Dyn),
    length(tstart_Dyn)
  )
)

# Create a pie chart
ggplot(count_data, aes(x = "", y = Count, fill = Variable)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_brewer(palette = "Set1") +
  theme(
    legend.text = element_text(size = 14)  # Adjust the size of legend text
  )




# Change colors of pie chart -- blue colors for increase, red for decrease 
blue_colors <- brewer.pal(9, "Blues")
red_colors <- brewer.pal(4, "Reds")

ggplot(count_data, aes(x = "", y = Count, fill = Variable)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  # scale_fill_manual(values = c("delta" = "#4292C6", 
  #                              "pT0" = "#DE2D26", 
  #                              "c" = "#6BAED6",
  #                              "beta" = "#377EB8", 
  #                              "tstart" = "#9ECAE1")) +
  scale_fill_manual(values = c("delta" = "#377EB8", 
                               "pT0" = "#DE2D26", 
                               "c" = "#4292C6",
                               "beta" = "#6BAED6", 
                               "tstart" = "#9ECAE1")) +
  theme(
    legend.text = element_text(size = 14)  # Adjust the size of legend text
  ) 



 #Check parameter values of best models 
B2AcReLogit50AllTotal <- read_csv("B2AcReLogit50AllTotal.csv")

BestID<- Summary_B2AcReLogit50Total$Run[1:rows]

BestModels<-NULL
for(i in 1:length(BestID)){
  x<-which(B2AcReLogit50AllTotal$Run==BestID[i])
  
  currentmodel<-B2AcReLogit50AllTotal[x, ]
  BestModels<-rbind(BestModels,currentmodel)
}

#delta
beta_delta_Dynamics_Rebound<-BestModels$estimate[which(BestModels$param=="beta_delta_Dynamics_Rebound")]


#log10pT0
beta_log10pT0_Dynamics_Rebound<-BestModels$estimate[which(BestModels$param=="beta_log10pT0_Dynamics_Rebound")]


#c
beta_c_Dynamics_Rebound<-BestModels$estimate[which(BestModels$param=="beta_c_Dynamics_Rebound")]

#log10beta
beta_log10beta_Dynamics_Rebound<-BestModels$estimate[which(BestModels$param=="beta_log10beta_Dynamics_Rebound")]

#tstart
beta_tstart_Dynamics_Rebound<-BestModels$estimate[which(BestModels$param=="beta_tstart_Dynamics_Rebound")]

#corr between c and beta 
corr_log10beta_c<-BestModels$estimate[which(BestModels$param=="corr_log10beta_c")]













