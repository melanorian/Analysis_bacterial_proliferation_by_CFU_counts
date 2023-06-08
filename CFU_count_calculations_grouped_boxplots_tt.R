#title: "CFU_count_calculations_grouped_boxplots.R"
#author: "Mela3ie Mendel"
#date: "23-3-2023"
#output: html_document
#update 8-6-23

# 1. Basic set up & loading data ----
# load libraries
library(ggplot2)
library(dplyr)
library(multcompView)
library(agricolae)
library(ggsignif)
library(ggpubr)
library("report")

# clean environment
rm(list=ls())

# Set workingdirectory
setwd("C:/Users/Mende012/Documents/Bioinformatics/CFU_counts/Analysis_bacterial_proliferation_by_CFU_counts")

# generate prefix string including date/ Initials
current_date <- gsub("-", "", as.character(Sys.Date()))# Get the current date as character
pre <- paste("MM", current_date, sep = "")

# filter for samples of interest 
pst_raw <- read.xlsx ("CFU_count_calculations_Pseudomonas_strains.xlsx", detectDates = T)

# 2. Select/ filter data ----
# 2.1 Pseudomonas strains as factors
pst_raw$Pseudomonas.strain <- as.factor(pst_raw$Pseudomonas.strain)
pst_raw$Date <- as.factor(pst_raw$Date)

# remove rows with NA (if extra rows are loaded)
# pst_raw <- pst_raw[1:210,]

# 2.2 log10 transformation of data for better presentation
MIN <- min(pst_raw$`Calculate.CFU/ml.per.1.cm2.of.harvested.leaf`[pst_raw$`Calculate.CFU/ml.per.1.cm2.of.harvested.leaf`!=0]) # set lowest CFU counts as 0
pst_raw$CFU.log10 <- log10(pst_raw$`Calculate.CFU/ml.per.1.cm2.of.harvested.leaf`+MIN)

# 2.3 make df containing only data I want to plot: Pst strain, CFU log10 transormed
Pst_strain <- as.factor(pst_raw$Pseudomonas.strain)
Pst_CFU <- pst_raw$CFU.log10
Pst_exp_date <- as.factor(pst_raw$Date)
Pst_OD <- pst_raw$OD6000.infiltrated.0dpi
Pst_comment <- pst_raw$Comment

simple_df <- data.frame(Pst_strain, Pst_OD, Pst_exp_date, Pst_CFU, Pst_comment)

# 2.4 Replace NAs in Pst_comment with "infiltrated"
simple_df$Pst_comment[is.na(simple_df$Pst_comment)] <- "infiltrated"

# 2.5 filtre for infiltrated samples
simple_df_un <- simple_df[simple_df$Pst_comment == "uninfiltrated", ]
simple_df <- simple_df[simple_df$Pst_comment == "infiltrated", ]
simple_df <- simple_df[!(simple_df$Pst_OD == "c"), ]

# 3. Basic Box Plots comparing Pst strains ----
# 3.1 subset into different OD
df_OD2 <- simple_df[(simple_df$Pst_OD == "0.2"),]
df_OD4 <- simple_df[(simple_df$Pst_OD == "0.4"),]
df_OD6 <- simple_df[(simple_df$Pst_OD == "0.6"),]


# 3.2 basic boxplot OD 0.2-0.6 - UPDATE OD in script if needed
 OD2 <- ggplot(df_OD2, aes(x=as.factor(Pst_strain), y=Pst_CFU, color = Pst_exp_date)) +
  geom_boxplot(#fill= c("seagreen4", "red")                # box colour
               outlier.colour = NA,           # Outliers color, 
               alpha=0.9)+                       # Box color transparency
  labs(x= "Pseudomonas strain", y = "log10 (CFU/cm2)") +
  
  geom_point(position = position_jitterdodge(0.1)) +
  
  guides(x = guide_axis(angle = 45 ))  +          # axis label
  theme_classic() +
  ylim(4, 8.5)

OD2


# basic boxplot OD 0.4
 OD4 <- ggplot(df_OD4, aes(x=as.factor(Pst_strain), y=Pst_CFU, color = Pst_exp_date)) +
  geom_boxplot(#fill= c("seagreen4", "red")                # box colour
               outlier.colour = NA,           # Outliers color, 
               alpha=0.9)+                       # Box color transparency
  labs(x= "Pseudomonas strain", y = "log10 (CFU/cm2)") +
  
  geom_point(position = position_jitterdodge(0.1)) +
  
  guides(x = guide_axis(angle = 45 ))  +          # axis label
  theme_classic() +
  ylim(4, 8.5)

OD4


# basic boxplot OD 0.6
 OD6 <- ggplot(df_OD6, aes(x=as.factor(Pst_strain), y=Pst_CFU, color = Pst_exp_date)) +
  geom_boxplot(#fill= c("seagreen4", "red")                # box colour
               outlier.colour = NA,           # Outliers color, 
               alpha=0.9)+                       # Box color transparency
  labs(x= "Pseudomonas strain", y = "log10 (CFU/cm2)") +
  
  geom_point(position = position_jitterdodge(0.1)) +
  
  guides(x = guide_axis(angle = 45 ))  +          # axis label
  theme_classic() +
  ylim(4, 8.5)

OD6

ggarrange(OD2, OD4, OD6, ncol = 3)

# 4. Boxplot grouped by strain strain over time (UPDATE: strains if needed)----

# 4.1 Get unique subsets based on Pst_strain
# simple_df <- simple_df[!(simple_df$Pst_OD == "c"),]

unique_strains <- unique(simple_df$Pst_strain)

# 4.2 Define the function to generate and save plots for each subset

co <- c("#44AA99", "#E69F00") # define 2 colours 
generate_plot <- function(df, strain) {
  subset_df <- df[df$Pst_strain == strain, ]
  
  g <- ggplot(subset_df, aes(x = as.factor(Pst_OD), y = Pst_CFU, color = Pst_exp_date)) +
    geom_boxplot(outlier.colour = "red", alpha = 0.9) +
    labs(x = strain, y = "log10 (CFU/cm2)") +
    geom_point(position = position_jitterdodge(0.1)) +
    theme_classic() +
    guides(x = guide_axis(angle = 45)) +
    scale_fill_manual(values = co) +
    scale_color_manual(values = co) +
    ylim(4, 9)
  
  ggsave(paste0(pre, "boxplot_", strain, ".svg"), plot = g, width = 3, height = 4, units = "in", dpi = 300)
}

# 4.3 Call the function using lapply for each unique strain
lapply(unique_strains, generate_plot, df = simple_df)

# Construction site start _______________________________________________________________

# 5. Statistical analysis that compares for each OD, Pst strain 0dpi and 2dpi (or days of interest)
# 5. Statistical analysis using a two-sided t-test ---- 

# 5.1 DEFINE VARIABLES OF INTEREST
OD <- unique(simple_df$Pst_OD) # extract the ODs from the orginal df
d1 <- "2023-01-19"  # Defines the reference dpi (usally 0dpi), epxeriment_date column as input  
d2 <- "2023-01-21"  # Defines dpi of interest, usually 2dpi, epxeriment_date column as input

# 5.2 Define the tt function: It takes as input the df of 1 bacterial strain, a specific OD
#     and two dates that are to be compared

tt <- function(df_Pst, OD, date1, date2) {
  a <- as.data.frame(df_Pst[(df_Pst$Pst_OD == OD), ])
  gr1 <- a[(a$Pst_exp_date == date1), ]
  gr2 <- a[(a$Pst_exp_date == date2), ]
  b <- t.test(gr1$Pst_CFU, gr2$Pst_CFU, alternative = "two.sided", var.equal = TRUE)
  return(b$p.value)
  }
  
# 5.3 Create a list of subsets using lapply and split
subset_list <- lapply(unique_strains, function(strain) {
    subset_df <- simple_df[simple_df$Pst_strain == strain, ]
    split(subset_df, subset_df$Pst_OD)
  })  

# 5.4 Perform t-tests on each subset using lapply
result_list <- lapply(subset_list, function(subset) {
    lapply(subset, function(df) {
      pval <- tt(df_Pst = df, OD = df$Pst_OD[1], date1 = d1, date2 = d2)
      data.frame(strain = df$Pst_strain[1], dpi = df$Pst_OD[1], p_value = pval)
    })
  })

# 5.6 Combine the results into a single dataframe
df_t_res <- do.call(rbind, unlist(result_list, recursive = FALSE)) # recursive = F: nested list, only first level is flattened into vector (1 df/ tt)

# 5.7 Create an empty significance column
df_t_res$significance <- ""
  
# 5.8 Update the significance column based on the p-value
df_t_res[df_t_res$p_value > 0.05, "significance"] <- "n.s"
df_t_res[df_t_res$p_value < 0.05, "significance"] <- "*"
df_t_res[df_t_res$p_value < 0.01, "significance"] <- "**"
df_t_res[df_t_res$p_value < 0.001, "significance"] <- "***"
  
# 5.9 safe the output as a .txt file
sink(paste(pre, sep = "","_tt_summary_CUF_Pst_strains.txt"))
print(df_t_res)
sink()
