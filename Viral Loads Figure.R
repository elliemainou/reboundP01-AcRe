# Histogram of c for B2AcReLogit50

library(deSolve)
library(ggplot2)
library(readr)
library(stringr) 
library(patchwork)
library(tidyverse)
setwd("/Users/ejm6007/Desktop/NLME/B2AcReLogit50")


PVL_B2_02172020_CG <- read_csv("~/PVL_B2_02172020_CG.csv")

ggplot(PVL_B2_02172020_CG, aes(x = week_infection, y = log10vl, color = factor(animal_id))) +
  geom_line() +
  labs(x = "Week of infection", y = "log10 viral load", color = "Animal ID") +
  theme_minimal() +
  theme(axis.text = element_text(size = 16, family = "Times New Roman"),  # Adjust the size and font for axis text
        axis.title = element_text(size = 18, family = "Times New Roman"),  # Set the font for axes titles
        plot.title = element_text(size = 18, face = "bold", family = "Times New Roman"), 
        legend.text = element_text(size = 16),  # Adjust legend text size
        #legend.title = element_text(size = 18), 
        text = element_text(family = "Times New Roman"))


ggplot(PVL_B2_02172020_CG, aes(x = week_infection, y = log10vl, color = factor(animal_id))) +
  geom_point() +
  geom_line() +
  labs(x = "Weeks post infection", y = "Log10 viral load(RNA copies/mL)",  color = "Animal ID") +
  theme_minimal() +
  theme(axis.text = element_text(size = 20, family = "Times New Roman"),  # Adjust the size and font for axis text
        axis.title = element_text(size = 22, family = "Times New Roman"),  # Set the font for axes titles
        plot.title = element_text(size = 18, face = "bold", family = "Times New Roman"), 
        legend.text = element_text(size = 16),  # Adjust legend text size
        #legend.title = element_text(size = 18), 
        text = element_text(family = "Times New Roman"))
