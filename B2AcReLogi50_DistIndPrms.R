# Histogram of c for B2AcReLogit50

library(deSolve)
library(ggplot2)
library(readr)
library(stringr) 
library(patchwork)
library(tidyverse)
setwd("/Users/ejm6007/Desktop/NLME/B2AcReLogit50")

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


par(mfrow=c(3,3))

## Acute= red
## Rebound=blue

# Install and load necessary packages if not already installed
# install.packages(c("patchwork", "tidyverse"))


# Create an empty list to store the individual plots
plots_list <- list()

# Loop through each iteration
for (i in 1:length(Run)){
  current <- BestModels[which(BestModels$x == Run[i]), ]
  Acute <- current[which(current$Dynamics == "Acute"), ]
  colnames(Acute)[4:8] <- c("deltaA", "log10pT0A", "cA", "log10betaA", "tstartA")
  Rebound <- current[which(current$Dynamics == "Rebound"), ]
  colnames(Rebound)[4:8] <- c("deltaR", "log10pT0R", "cR", "log10betaR", "tstartR")
  
  # Create individual plots and store in the list
  plot <- ggplot() +
    geom_histogram(data = gather(Acute, key = "group", value = "value", deltaA), 
                   aes(x = value, fill = "Acute"),
                   position = "identity", alpha = 0.5, bins = 10) +
    geom_histogram(data = gather(Rebound, key = "group", value = "value", deltaR), 
                   aes(x = value, fill = "Rebound"),
                   position = "identity", alpha = 0.5, bins = 10) +
    labs(title = paste("Model",  Run[i]),
         x = "δ",
         y = "Frequency") +
    theme_minimal() +
    scale_fill_manual(values = c("red", "blue")) +  # Removed the name argument
    theme(axis.text = element_text(size = 18, family = "Times New Roman"),
          axis.title = element_text(size = 20, family = "Times New Roman"),
          plot.title = element_text(size = 24, face = "bold", family = "Times New Roman"), 
          legend.position = "none",  # Remove the legend
          text = element_text(family = "Times New Roman"))
  
  plots_list[[i]] <- plot
}

# Arrange the plots into a 4x3 grid using patchwork
combined_plot <- wrap_plots(plots_list, nrow = 4)

# Print or save the combined plot
print(combined_plot)


# Load ggplot2 package
library(ggplot2)

# Create ggplot with scatterplot and faceting by 'x'
ggplot(data = BestModels, aes(x = delta, y = log10pT0)) +
  geom_point() +  # Create scatterplot with points
  facet_wrap(~ x) +  # Create separate panel for each unique value of 'x'
  labs(
    x = "Delta",  # Label for the x-axis
    y = "log10(pT0)",  # Label for the y-axis
    title = "Scatterplot of Delta vs. log10(pT0) by x"  # Title of the plot
  ) +
  theme_minimal()  # Optional: a cleaner theme for visualization

# Load ggplot2 and dplyr packages

# Calculate correlations for each group of 'x'
correlations <- BestModels %>%
  group_by(x) %>%
  summarize(correlation = cor(delta, log10pT0))

# Create ggplot with scatterplot, faceting by 'x', and color by 'Dynamics'
scatter_plot <- ggplot(data = BestModels, aes(x = delta, y = log10pT0, color = Dynamics)) +
  geom_point() +  # Create scatterplot with points
  facet_wrap(~ x) +  # Separate panel for each unique value of 'x'
  labs(
    x = "Delta",
    y = "log10(pT0)",
    title = "Scatterplot of Delta vs. log10(pT0) by x",
    color = "Dynamics"  # Label for the color legend
  ) +
  theme_minimal()  # Use a cleaner theme for visualization

# Modify the color scheme to have "Acute" as red and others default
scatter_plot <- scatter_plot +
  scale_color_manual(values = c("Acute" = "red", "Rebound" = "blue"))  # Define custom colors

# Add the calculated correlations to the plot as annotations
scatter_plot_with_cor <- scatter_plot +
  geom_text(
    data = correlations,
    aes(x = Inf, y = Inf, label = sprintf("r = %.2f", correlation)),
    hjust = 1.1,  # Horizontal adjustment to position the text
    vjust = 1.1,  # Vertical adjustment to position the text
    size = 3.5,  # Size of the text
    color = "black"  # Color of the text (different from scatter points for contrast)
  )

# Display the final plot with faceting, color by 'Dynamics', and correlations
print(scatter_plot_with_cor)


#####################################
# Subset the data to include only points where 'Dynamics' is 'Acute'
acute_data <- BestModels %>%
  filter(Dynamics == "Rebound")

# Calculate correlations for 'Acute' data by 'x'
acute_correlations <- acute_data %>%
  group_by(x) %>%
  summarize(correlation = cor(delta, log10pT0))

# Create a scatterplot with 'delta' and 'log10pT0', faceting by 'x'
scatter_plot <- ggplot(data = acute_data, aes(x = delta, y = log10pT0)) +
  geom_point(color = "red") +  # Points are red since all are 'Acute'
  facet_wrap(~ x) +  # Create separate panels for each unique 'x'
  labs(
    x = "Delta",
    y = "log10(pT0)",
    title = "Scatterplot of Delta vs. log10(pT0) by x (Acute Only)"
  ) +
  theme_minimal()  # Minimal theme for cleaner visualization

# Add the calculated correlations to the plot as annotations
scatter_plot_with_cor <- scatter_plot +
  geom_text(
    data = acute_correlations,
    aes(x = Inf, y = Inf, label = sprintf("r = %.2f", correlation)),
    hjust = 1.1,  # Horizontal adjustment to position the text
    vjust = 1.1,  # Vertical adjustment to position the text
    size = 3.5,  # Size of the text
    color = "black"  # Color for annotation text
  )

# Display the plot with correlations
print(scatter_plot_with_cor)
