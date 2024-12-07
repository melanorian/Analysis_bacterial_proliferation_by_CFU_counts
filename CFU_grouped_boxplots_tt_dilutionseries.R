#title: "CFU_count_calculations_grouped_boxplots.R"
#author: "Mela3ie Mendel"
#date: "23-3-2023"
#output: html_document
#update 09-10-24

# clean environment
rm(list=ls())

# 1. Basic set up & loading data ----
# load libraries
library(ggplot2)
library(dplyr)
library(multcompView)
library("agricolae")
library(ggsignif)
library("ggpubr")
library("report")
library("readxl")

# Set workingdirectory
variance_for_tt <- TRUE
SUF <- "_rm0_nomin_nooutlier"
setwd("/home/melanie/working_directory/D36E_assays/CFU_data")
dir_out <- "/home/melanie/working_directory/D36E_assays/CFU_output/20240920/t0_t2_t6" 

# generate prefix string including date/ Initials
current_date <- gsub("-", "", as.character(Sys.Date()))# Get the current date as character
pre <- paste("MM", current_date, sep = "")

# filter for samples of interest 
pst_raw <- read_xlsx("2025092024_CFU_count_calculations_dilutions.xlsx")

# 2. Select/ filter data ----
# 2.1 Pseudomonas strains as factors
pst_raw$Pseudomonas.strain <- as.factor(pst_raw$`Pseudomonas strain`)
pst_raw$Date <- as.factor(pst_raw$dpi)

# remove rows with NA (if extra rows are loaded)
# pst_raw <- pst_raw[1:210,]
# rm if 0
zeros <- pst_raw[!pst_raw$`Calculate CFU/ml per 1 cm2 of harvested leaf` == 0 ,]
cat("nr removed rows due to zero-value",nrow(pst_raw) - nrow(zeros))
pst_raw <- zeros
rm(zeros)

# 2.2 log10 transformation of data for better presentation
#MIN <- min(pst_raw$`Calculate CFU/ml per 1 cm2 of harvested leaf`[pst_raw$`Calculate CFU/ml per 1 cm2 of harvested leaf`!=0]) # set lowest CFU counts as 0
pst_raw <- pst_raw[!is.na(as.numeric(pst_raw$`Calculate CFU/ml per 1 cm2 of harvested leaf`)), ] # remove rows with na
pst_raw$CFU.log10 <- log10(as.numeric(pst_raw$`Calculate CFU/ml per 1 cm2 of harvested leaf`)) #+MIN)

# 2.3 make df containing only data I want to plot: Pst strain, CFU log10 transormed
Pst_strain <- as.factor(pst_raw$Pseudomonas.strain)
Pst_CFU <- pst_raw$CFU.log10
Pst_exp_date <- as.factor(pst_raw$Date)
Pst_OD <- pst_raw$`OD6000 infiltrated 0dpi`
Pst_comment <- pst_raw$Notes

simple_df <- data.frame(Pst_strain, Pst_OD, Pst_exp_date, Pst_CFU)

# 2.5 filtre for infiltrated samples
simple_df_un <- simple_df[simple_df$Pst_exp_date == "i", ] # infiltrate samples
simple_df <- simple_df[!(simple_df$Pst_exp_date == "i"), ] # in planta samples

# Remove 6 dpi samples, "bad" samples 
#simple_df <- simple_df[!(simple_df$Pst_exp_date == "6"), ]
#simple_df <- simple_df[!(simple_df$Pst_comment == "bad"), ]

df <- simple_df
# Remove outliers
remove_outliers <- function(df) {
df %>%
group_by(Pst_strain, Pst_OD, Pst_exp_date) %>%

filter(Pst_CFU >= (quantile(Pst_CFU, 0.25) - 1.5 * IQR(Pst_CFU)) &
 Pst_CFU <= (quantile(Pst_CFU, 0.75) + 1.5 * IQR(Pst_CFU)))
}

# Apply the function
cleaned_df <- remove_outliers(simple_df)
cat("nr removed samples:", nrow(simple_df) - nrow(cleaned_df))
simple_df <- cleaned_df

# 3. Basic Box Plots comparing Pst strains ----
# 3.2 basic boxplot OD 0.2-0.6 - UPDATE OD in script if needed
unique_OD <- unique(simple_df$Pst_OD)

# 4.2 Define the function to generate and save plots for each subset

generate_OD_plot <- function(df, OD) {
subset_df <- df[df$Pst_OD == OD, ]

g <- ggplot(subset_df, aes(x = as.factor(Pst_strain), y = Pst_CFU, color = Pst_exp_date)) +
geom_boxplot(outlier.colour = "red", alpha = 0.9) +
labs(x = OD, y = "log10 (CFU/cm2)") +
ggtitle(paste0("OD: ", OD)) +
geom_point(position = position_jitterdodge(0.1)) +
theme_classic() +
guides(x = guide_axis(angle = 45)) +
scale_fill_manual(values = co) +
scale_color_manual(values = co) +
ylim(2.5, 9)

return(g)
}

# Generate the plots using lapply
co <- c("#44AA99", "#E69F00", "#999999")

#simple_df <- simple_df[(simple_df$Pst_strain == "D36" | simple_df$Pst_strain == "DC3000"),]
OD_plot_list <- lapply(unique_OD, generate_OD_plot, df = simple_df)
OD_grid <- ggarrange(plotlist = OD_plot_list, ncol = length(unique_OD), nrow = 1)

# # Save the plots using purrr::map2 and a lambda function
# map2(.x = OD_plot_list, 
# .y = unique_OD, 
# .f = ~ggsave(filename = paste0(dir_out, pre, "_boxplot_", SUF, .y, ".svg"),
#       plot = .x, 
#       width = 3.5, 
#       height = 4, 
#       units = "in", 
#       dpi = 300)
# )

ggsave(filename = paste0(dir_out, pre, "_Lines_grid", SUF, ".svg"),
plot = OD_grid, 
width = 12, 
height = 4, 
units = "in", 
dpi = 300)

# 4. Boxplot grouped by strain strain over time (UPDATE: strains if needed)----

# 4.1 Get unique subsets based on Pst_strain
# simple_df <- simple_df[!(simple_df$Pst_OD == "c"),]

unique_strains <- unique(simple_df$Pst_strain)

# 4.2 Define the function to generate and save plots for each subset

generate_strain_plot <- function(df, strain) {
subset_df <- df[df$Pst_strain == strain, ]

g <- ggplot(subset_df, aes(x = as.factor(Pst_OD), y = Pst_CFU, color = Pst_exp_date)) +
geom_boxplot(outlier.colour = "red", alpha = 0.9) +
labs(x = strain, y = "log10 (CFU/cm2)") +
ggtitle("Strain:", strain) +
geom_point(position = position_jitterdodge(0.1)) +
theme_classic() +
guides(x = guide_axis(angle = 45)) +
scale_fill_manual(values = co) +
scale_color_manual(values = co) +
ylim(2.5, 9)

return(g)
}

# Generate the plots using lapply
strain_plot_list <- lapply(unique_strains, generate_strain_plot, df = simple_df)
strain_grid <- ggarrange(plotlist = strain_plot_list, ncol = length(unique_strains), nrow = 1)


#################

generate_line_plot <- function(df, strain) {
  subset_df <- df[df$Pst_strain == strain, ]
  
  g <- ggplot(subset_df, aes(x = as.factor(Pst_OD), y = Pst_CFU, color = Pst_exp_date)) +
    geom_boxplot(outlier.colour = "red", alpha = 0.9) +
    labs(x = strain, y = "log10 (CFU/cm2)") +
    ggtitle("Strain:", strain) +
    geom_point(position = position_jitterdodge(0.1)) +
    theme_classic() +
  guides(x = guide_axis(angle = 45)) +
    scale_fill_manual(values = co) +
    scale_color_manual(values = co) +
    ylim(2.5, 9)
  
  return(g)
}

# Generate the plots using lapply
strain_plot_list <- lapply(unique_strains, generate_strain_plot, df = simple_df)
strain_grid <- ggarrange(plotlist = strain_plot_list, ncol = length(unique_strains), nrow = 1)


# # Save the plots using purrr::map2 and a lambda function
# map2(.x = strain_plot_list, 
# .y = unique_strains, 
# .f = ~ggsave(filename = paste0(dir_out, pre, "_boxplot_", SUF, .y, ".svg"),
#       plot = .x, 
#       width = 3, 
#       height = 4, 
#       units = "in", 
#       dpi = 300)
# )

ggsave(filename = paste0(dir_out, pre, "_strain_grid", SUF, ".svg"),
plot = strain_grid, 
width = 12, 
height = 4, 
units = "in", 
dpi = 300)

# 5. Statistical analysis that compares for each OD, Pst strain 0dpi and 2dpi (or days of interest)
# 5. Statistical analysis using a two-sided t-test ---- 

# 5.1 DEFINE VARIABLES OF INTEREST
OD <- unique(simple_df$Pst_OD) # extract the ODs from the orginal df
d1 <- "0"  # Defines the reference dpi (usally 0dpi), epxeriment_date column as input  
d2 <- "3"  # Defines dpi of interest, usually 2dpi, epxeriment_date column as input

# 5.2 Define the t-test function
f_tt <- function(df_Pst, OD, date1, date2, variance_equal) {
gr1 <- df_Pst[df_Pst$Pst_exp_date == date1, ]
gr2 <- df_Pst[df_Pst$Pst_exp_date == date2, ]
b <- t.test(gr1$Pst_CFU, gr2$Pst_CFU, alternative = "two.sided", var.equal = variance_equal)
return(b$p.value)
}

f_perform_tt <- function(subset_list, date1, date2, variance_equal) {

# 1. Perform t-tests on each subset using lapply
result_list <- lapply(subset_list, function(subset) {


lapply(subset, function(df) {
pval <- f_tt(df_Pst = df,
 #OD = df$Pst_OD[1], # start from first element
 date1 = date1, 
 date2 = date2, 
 variance_equal = variance_equal
 )
data.frame(strain = df$Pst_strain[1], dpi = df$Pst_OD[1], p_value = pval)
})
})

# 2. Combine the results into a single dataframe
df_t_res <- do.call(rbind, unlist(result_list, recursive = FALSE))

# 3. Create a significance column based on the p-values
df_t_res$significance <- ""
df_t_res[df_t_res$p_value > 0.05, "significance"] <- "n.s"
df_t_res[df_t_res$p_value < 0.05, "significance"] <- "*"
df_t_res[df_t_res$p_value < 0.01, "significance"] <- "**"
df_t_res[df_t_res$p_value < 0.001, "significance"] <- "***"

# 4. Update column names and row names
colnames(df_t_res) <- c("strain", "OD600", "p-val", "significance")
rownames(df_t_res) <- 1:nrow(df_t_res)

# 5. Save the output to a .txt file
tt <- ifelse(variance_equal, "_equalvar", "_welch")
output_file <- paste0(dir_out, pre, "_tt_twosided_t0_vs_t2", SUF, tt, ".txt")
sink(output_file)
print(df_t_res)
sink()

# Return the result dataframe
return(cat(output_file))
}

unique_strains <- c("D36", "DC3000")
subset_list <- lapply(unique_strains, function(strain) {
  subset_df <- simple_df[simple_df$Pst_strain == strain, ]
  split(subset_df, list(subset_df$Pst_OD, subset_df$Pst_strain))
})

cleaned_subset_list <- lapply(subset_list, function(sublist) {
  Filter(function(df) nrow(df) > 0, sublist)
})
# # 5.3 Create a list of subsets using lapply and split
# subset_list <- lapply(unique_strains, function(strain) {
# subset_df <- simple_df[simple_df$Pst_strain == strain, ]
# split(subset_df, subset_df$Pst_OD)
# })

res_tt <- f_perform_tt(subset_list = cleaned_subset_list, 
date1 = d1, 
date2 = d2, 
variance_equal = variance_for_tt
)


### ANOVA + PstHoc for subsets ###### ----

# A) ANOVA 
subset_df <- simple_df[simple_df$Pst_OD == "0.4", ]

subset_df <- subset_df[subset_df$Pst_strain == "AvrE1",]
subset_df$ID <- paste(subset_df$Pst_strain, subset_df$Pst_exp_date, sep = "_")

anova_pst <- aov(Pst_CFU ~ ID, data = subset_df)

{summary(anova_pst)
  summary.lm(anova_pst)}

# B) PostHoc

# 2. Posthoc to see which groups are significantly different - generate labels

# posthoc
TUKEY <- TukeyHSD(anova_pst)

generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$ID = rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$ID) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset

LABELS <- generate_label_df(TUKEY , "ID")
yvalue <- aggregate(Pst_CFU ~ ID, data = subset_df, mean) # obtain letter position for y axis using means
final <- merge(yvalue,LABELS, by = "ID") 

output_file <- paste0(dir_out, "/", pre, "_TUKEY_t0_to_t6_OD04_AvrE1", ".txt")
sink(output_file)
print(final)
sink()

stat_labels <- as.vector(LABELS$Letters) # extract letter column as a vector to use for plot
