---
# title: "Final Analysis CFU count DC3000 Effector library"
# author: "Melanie Mendel"
# date: "07-01-2023"
# output: html_document
# update: "06-06-2023"
---
  
# 1. Basic settings & data import ----

# clean up
rm(list=ls())

# set wd for whole document
setwd("C:/Users/Mende012/Documents/Bioinformatics/CFU_counts/Analysis_bacterial_proliferation_by_CFU_counts")

# generate prefix string including date/ Initials
current_date <- gsub("-", "", as.character(Sys.Date()))# Get the current date as character
pre <- paste("MM", current_date, sep = "")

# load required packages
library(dplyr)
library(multcompView)
library(agricolae)
library(ggsignif)
library(ggpubr)
library("report")
library("ggpubr")
library(tidyverse)
library(plotly)
library(openxlsx)

# filter for samples of interest 
pst_raw <- read.xlsx ("CFU_count_calculations_DC3000_library.xlsx", detectDates = T)
head(pst_raw)

# 2. Generate workable df ----
# 2.1 filter out all samples at t2
pst_raw <- pst_raw %>% filter(dpi =="2")
pst_raw[is.na(pst_raw)] <- 0

# 2.2 Pseudomonas strains as factors
pst_raw$Pseudomonas.strain <- as.factor(pst_raw$Pseudomonas.strain)

# 2.3 log10 transformation of data for better presentation
MIN <- min(pst_raw$`Calculate.CFU/ml.per.1.cm2.of.harvested.leaf`[pst_raw$`Calculate.CFU/ml.per.1.cm2.of.harvested.leaf`!=0]) # set lowest CFU counts as 0
pst_raw$CFU.log10 <- log10(pst_raw$`Calculate.CFU/ml.per.1.cm2.of.harvested.leaf`+MIN) # log 10 transformation

# 2.4 make df containing only data I want to plot: Pst strain, CFU log10 transormed
Pst_strain <- pst_raw$Pseudomonas.strain
Pst_CFU <- pst_raw$CFU.log10
Date <- as.factor(pst_raw$Date) ###

# 2.5 make simple df
simple_df <- data.frame(Pst_strain, Pst_CFU, Date)
head(simple_df)



# ------------------------- Data not normalised/log10 transformed -----------------------
# 3. Box Plots comparing Pst strains, 2 dpi, no stats, save as pdf ----
# If printed: CHANGE FILE LOCATION


#pdf(file = "C:/Users/Mende012/Documents/Bioinformatics/CFU_counts/MM20230107_CFU_count_raw.pdf",
#     width = 8,
#     height = 5)

# 3.1 basic boxplot
g1 <-  ggplot(simple_df, aes(x=as.factor(Pst_strain), y=Pst_CFU)) +
  
  geom_boxplot(fill="white",                 # box colour
               outlier.colour = "white",     # Outliers color, 
               alpha=0.9) +                  # Box color transparency
  
  theme_bw() +  # make the bg white
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove background, frame
        axis.line = element_line(colour = "black")) +
  
  xlab("Pseudomonas strain") +                   # axis label
  ylab("log10 (CFU/cm2)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # turn the tick marks on the x axis 45 degree
  geom_jitter(shape=16, position=position_jitter(0.1),
              aes(colour = factor(Date))) +
  stat_boxplot(geom = "errorbar", # Error bars
                width = 0.2)

g1

# Run dev.off() to create the file!
dev.off()

# 4. Statistical analysis across samples using ANOVA+TUKEY ----
# 4.1. ANOVA to analyse whether there are differences between the groups
#      H0 All means are the same 
#      H1 at least one mean is different

anova_pst <- aov(Pst_CFU ~ Pst_strain, data = simple_df)

  {summary(anova_pst)
    summary.lm(anova_pst)}

# 4.2 TUKEY Posthoc to see which groups are significantly different 
TUKEY <- TukeyHSD(anova_pst)

# 4.3 Summarize TUKEY test results for saving & plotting

# 4.3.1 write a funktion to generate labels for the test results
generate_label_df <- function(TUKEY, variable){
    
    # Extract labels and factor levels from Tukey post-hoc results
    Tukey.levels <- TUKEY[[variable]][,4]
    Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
    
    # order labels according to boxplot:
    Tukey.labels$Pst_strain = rownames(Tukey.labels)
    Tukey.labels=Tukey.labels[order(Tukey.labels$Pst_strain) , ]
    return(Tukey.labels)
    }

# 4.3.2 Apply the function on the df
LABELS <- generate_label_df(TUKEY , "Pst_strain")
names(LABELS) <- c("Letters","Pst_strain")

# 4.3.3 safe TUKEY test results
sink(paste(pre, sep = "","CFU_count_DC3000_library_TUKEY_nonormal.txt"))
print(LABELS)
sink()
  
# 4.3.4 obtain letter position for y axis using means
yvalue <- aggregate(Pst_CFU ~ Pst_strain, data = simple_df, mean) 
LABELS$Pst_strain <- (unlist(lapply(strsplit(LABELS$Pst_strain,'[_]'),function(x){x[[1]]})))

# 4.3.5 generate df with results TUKEY test
final <- merge(yvalue,LABELS) 

# 4.3.6 For ploting extract label vector
stat_labels <- as.vector(LABELS$Letters)
  

# 5. Box Plots comparing Pst strains, 2 dpi, TUKEY results, save as pdf (default) ----

pdf(file = paste(pre, sep = "","_CFU_count_raw_TUKEY.pdf"),
     width = 10,
     height = 6)

g_2 <-  ggplot(simple_df, aes(x=as.factor(Pst_strain), y=Pst_CFU)) + 
  
  geom_boxplot(fill="white",                 # box colour
               outlier.colour = "white",     # Outliers color, 
               alpha=0.9) +                  # Box color transparency
  
  theme_bw() +  # make the bg white
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove background, frame
        axis.line = element_line(colour = "black")) +
  
  xlab("Pseudomonas strain") +                   # axis label
  ylab("log10 (CFU/cm2)") +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # turn the tick marks on the x axis 45 degree +
  expand_limits( y = c(0, 10)) +
  
  geom_jitter(shape=16, position=position_jitter(0.1),
              aes(colour = factor(Date))) +
  stat_boxplot(geom = "errorbar", # Error bars
                width = 0.2) +
  geom_text(data = final, aes(y=Pst_CFU, label = Letters), 
                              vjust= -7 ,hjust= 0.5)


g_2

# Run dev.off() to create the file!
dev.off()



# ------------------------- Data normalization/log10 transformed  -----------------------

# 6. normalise data by corresponding negative control (CHOOSE MEAN/MEDIAN) ----
#    here: D36E
# 6.1 extract the average CFU of the negative control for normalization df
# 6.1.1 List all the unique dates and store as a vector
date <- unique(simple_df$Date)
date <- as.vector(unlist(date))

# 6.1.2 write function to extract dates and calculate the mean

n_control_mean <- function(df, Pst_control_strain, date) {
  
 D36E <- subset(df, Pst_strain == Pst_control_strain)
 D36E_date <- subset(D36E, date == Date)
 D36E_mean_date <- mean(D36E_date$Pst_CFU) # CHANGE POSSIBLE MEAN/MEDIAN
 
 return(D36E_mean_date)
 
}

# 6.1.3 apply the function and extract mean of D36E control for the separate 
#       dates and store as a vector

n_control_mean <- lapply(date, n_control_mean, df = simple_df, Pst_control_strain = "D36")
n_control_mean <- as.vector(unlist(n_control_mean))
unlist(n_control_mean)

# 6.1.4 create df from control means
df_control_mean <- data.frame(date, n_control_mean)


# 6.2 normalise CFU data points using the mean/median of the control
# 6.2.1 write function to normalise CFU by date

general_normalise <- function(CFU_count, CFU_count_control){
  norm <- CFU_count - CFU_count_control
  return(norm)
}

# 6.2.2 Write function to extract values by date, and normalise them usign 6.2.1

normalise <- function(df, date){
  
  CFU_date <- subset(df, Date == date)
  CFU_date_Pst_CFU <- CFU_date$Pst_CFU
  
  date_mean <- df_control_mean[df_control_mean$date == date ,]
  date_mean <- as.numeric(date_mean$n_control_mean)
  
  date_norm <- as.vector(general_normalise(CFU_count = CFU_date_Pst_CFU, CFU_count_control = date_mean))
  }


# 6.2.3 Apply nested normalise function for all the dates of experiments
#       unlist and save as vector
CFU_norm <- as.vector(unlist(lapply(date, normalise, df = simple_df)))

# 6.3 Safe normalized data in a new simple_df
simple_df$Pst_CFU_norm <- CFU_norm

# 6.4 Create simple_df_norm that only includes the normalised CFU counts for plot
# 6.4.1 Extract individual vectors
Pst_strain_n <- simple_df$Pst_strain
Pst_CFU_n <- simple_df$Pst_CFU_norm
Date <- as.factor(simple_df$Date)

# 6.4.2 summarise the vectors in a df
simple_df_norm <- data.frame(Pst_strain_n, Pst_CFU_n, Date)
head(simple_df_norm)

# Optional: If sampples should be removed:
# e.g. remove D36E control after normalisation
# simple_df_norm <- subset(simple_df_norm, Pst_strain != "D36")

# 7. Statistical analysis after normalisation samples using ANOVA+TUKEY ----
# 7.1 ANOVA to analyse whether there are differences between the groups
#     (Details see 4.) 

anova_pst_norm <- aov(Pst_CFU_n ~ Pst_strain_n, data = simple_df_norm)

  {summary(anova_pst_norm)
    summary.lm(anova_pst_norm)}

# 7.2 TUKEY test 
TUKEY_n <- TukeyHSD(anova_pst_norm)
  
  generate_label_df <- function(TUKEY, variable){
    
    # Extract labels and factor levels from Tukey post-hoc
    Tukey.levels <- TUKEY[[variable]][,4]
    Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
    
    #I need to put the labels in the same order as in the boxplot :
    Tukey.labels$Pst_strain = rownames(Tukey.labels)
    Tukey.labels=Tukey.labels[order(Tukey.labels$Pst_strain) , ]
    return(Tukey.labels)
  }
  
# Apply the function on my dataset

  LABELS_n <- generate_label_df(TUKEY_n , "Pst_strain_n")
  names(LABELS_n) <- c("Letters_n","Pst_strain_n")
  yvalue_n <- aggregate(Pst_CFU_n ~ Pst_strain_n, data = simple_df_norm, mean) # obtain letter position for y axis using means
  LABELS_n$Pst_strain_n <- (unlist(lapply(strsplit(LABELS_n$Pst_strain_n,'[_]'),function(x){x[[1]]})))
  final_n <- merge(yvalue_n,LABELS_n) 
  final_n$Pst_CFU_n <- final_n$Pst_CFU_n + 1 # add +1 so that the label is not on the plot (e.g. outliers)
  stat_labels_n <- as.vector(LABELS_n$Letters_n) # extract letter column as a vector to use for plot

  # 4.3.3 safe TUKEY test results
  sink(paste(pre, sep = "","CFU_count_DC3000_library_TUKEY_normalised.txt"))
  print(LABELS_n)
  sink()
 
# 8. Box Plots comparing Pst strains, 2 dpi, TUKEY results, normalised save as pdf (default) ----
  
  pdf(file = paste(pre, sep = "","_CFU_count_TUKEY_normalised.pdf"),
      width = 10,
      height = 6)
  
  g_3 <-  ggplot(simple_df_norm, aes(x=as.factor(Pst_strain_n), y=Pst_CFU_n)) + 
    
    geom_boxplot(fill="white",                 # box colour
                 outlier.colour = "white",     # Outliers color, 
                 alpha=0.9) +                  # Box color transparency
    
    theme_bw() +  # make the bg white
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.border = element_blank(), # remove background, frame
          axis.line = element_line(colour = "black")) +
    
    xlab("Pseudomonas strain") +                   # axis label
    ylab("log10 (CFU/cm2)") +
    
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # turn the tick marks on the x axis 45 degree +
    expand_limits( y = c(0, 5)) +
    
    geom_jitter(shape=16, position=position_jitter(0.1),
                aes(colour = factor(Date))) +
    stat_boxplot(geom = "errorbar", # Error bars
                 width = 0.2) +
    geom_text(data = final_n, aes(y=Pst_CFU_n, label = Letters_n), 
              vjust= -7 ,hjust= 0.5)
  
  
  
  g_3
  
# Run dev.off() to create the file!
  dev.off()
  
  
  
  

# alternative stats Dunnett's test
```{r}
library("DescTools")
dt_no_anova <- DunnettTest(Pst_CFU_n ~ Pst_strain_n, data = simple_df_norm, control = "D36", conf.level = 0.95)

```


#DunnettTest(log2Norm_Fc ~ Pst_strain_n, data = FC[FC$Gene==gene & FC$Time==time,],
 #           control = "D36", conf.level = 0.95)"

# plotting after normalisation

```{r}
# add objects to overlay horizontal lines on plot
D36E_df <- simple_df_norm[(simple_df_norm$Pst_strain_n == "D36"),]
D36E_median <- median(D36E_df$Pst_CFU_n)
sd_plus <- D36E_median + sd(D36E_df$Pst_CFU_n)
sd_minus <- D36E_median - sd(D36E_df$Pst_CFU_n)

cbp1 <- c("#999999", "#44AA99", "#E69F00", "#661100", "#56B4E9", "#000000", "#117733", "#332288", "#F0E442", "#0072B2", "#D55E00", "#999933", "#CC79A7", "#E69F00")

# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# create pdf
pdf(file = "C:/Users/Mende012/Documents/Bioinformatics/CFU_counts/MM20230107_CFU_count_norm_D36E.pdf",
     width = 8,
     height = 5)

# basic boxplot
g2 <-  ggplot(
  
  #input data
  simple_df_norm, aes(x=as.factor(Pst_strain_n), y=Pst_CFU_n)) +
  
  # generate basic boxplot
  geom_boxplot(fill="white",                 # box colour
               outlier.colour = "white",     # Outliers color, 
               alpha=0) +                    # Box color transparency
  
  # overlay with jitter
  geom_jitter(shape=16, position=position_jitter(0.1),
              aes(colour = factor(Date))) +
  scale_colour_manual(values = cbp1) +
  
  # define the theme of the boxplot
  theme_bw() +  # make the bg white
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove background, frame
        axis.line = element_line(colour = "black")) +
  
  # label the axises 
  xlab("Pseudomonas strain") +                
  ylab("log10 fold change (CFU/cm2)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # turn the tick marks on the x axis 45 degree 
  

  # define axis limits if needed
  expand_limits( y = c(0, 5)) +
  
  # add statistical information: error bars, statistical results
  stat_boxplot(geom = "errorbar", # Error bars
                width = 0.2)  +
  
  # ANOVA/TUKEY test labels
  #geom_text(data = final_n, aes(y=Pst_CFU_n, label = Letters_n), 
                              #vjust= -10 ,hjust= 0.5) 
  
  # significance for single comparisons
  #geom_signif(comparisons = list(c("D36", "D36+AvrE1"),
   #                              c("D36", "D36+HopM1"),
    #                             c("D36", "D36+HopU1"),
     #                            c("D36", "DC3000")),
      #        map_signif_level = TRUE,
       #       y_position = c(3, 3.5, 4, 4.5))+
  
              
    ggsignif::geom_signif(
    comparisons = list(c("D36", "D36+AvrE1"),
                       c("D36", "D36+HopM1"),
                       c("D36", "D36+HopU1"),
                       c("D36", "DC3000")),
    annotations = c("***", "***", "**", " ***"),
    y_position = c(3, 3.5, 4, 4.5) 
  
  )
    
 # add ribbon to indicate st error of D36E
 # geom_ribbon(aes(ymin = sd_minus, ymax = sd_plus), alpha = 1) +  

# close pdf
ggsave("MM20230107_CFU_RLU_simple.svg",width = 7, height = 4)
g2
dev.off()

g2

# for adding dunnets test results see: https://jmbuhr.de/dataintro/freestyle.html
```

# for summary heat map extract medians of all samples
```{r}

# list all unique strains
strains <- unlist(unique(simple_df_norm$Pst_strain_n))

# subset normalized simple_df by single vectors and calculate mean/sample

#strain <- "D36+HopA1"
#rm(ls=strain)
median_values <- lapply(strains, function(x){
data <- simple_df_norm[simple_df_norm$Pst_strain_n == x,]
medians <- median(data$Pst_CFU_n)
return(medians)
})

# generate small df to store information
strains <- as.character(strains)
df_median <- cbind("Pst_strain" = strains, "median" = median_values)

# export 
write.csv(df_median, "./Median_CFU_counts.csv")
```

# -------------------- Data normalization/log10 transformed/scale by DC3000 --------------
# 9. Scale by positive control DC3000 ----
# 9.1 write function to show  ratio between median/median positive control
#     and sample measurement, more explanation in 7.

## write a function to identify the median/average of the positive control
# identify unique date and sort them to ensure order
date <- unique(simple_df_norm$Date)
date <- as.vector(sort(date))

# apply the function over all the assay s 
contr_median <-   lapply(date, function(x){
  a <- as.data.frame(simple_df_norm[simple_df_norm$Date == x,]) 
  b <- as.data.frame(a[a$Pst_strain_n == "DC3000",])
  contr_median <- median(b$Pst_CFU_n)
  return(contr_median)
})     

## make small df to store median/median results per assay
df_ctrmed <- data.frame(date, unlist(contr_median))
colnames(df_ctrmed) <- c("Date", "median_positive_control") #change median/median
rm(ls = contr_median)

# use the function with lapply to repeat the same action for different dates
rows <-   lapply(date, function(x){
  
  rows <- nrow(simple_df_norm[simple_df_norm$Date == x,])
  a <- df_ctrmed[df_ctrmed$Date == x,]
  re_assay_date <- unlist(rep(a$Date, times = rows))
  re_median <- unlist(rep(a$median_positive_control, times = rows))
  b <- data.frame(unlist(re_assay_date), unlist(re_median))
  return(b)
  }
  ) %>% bind_rows()

# 9.2 reorganise and merge df
# rename columns in a medianingful way 
colnames(rows) <- c("Date", "median_positive_control") #change name median/median

# order df so that it is also ordered by plateID 
df <- simple_df_norm[order(simple_df_norm$Date, decreasing = FALSE),]

# add additional column to original df
df["median_positive_control"] <- rows$median_positive_control

# 9.3 clean up 
rm(ls=df_ctrmed)
rm(ls=rows)


# 9.4 use the median/median to scale the measurements
# 9.4.1 write a function to calculate how the measurements compare to the positive control 
# calculate this as a % of the positive_control_median

prc_pos_avg <- (df$Pst_CFU_n/df$median_positive_control)*100
df["prc_pos_avg"] <- prc_pos_avg


# 10. Statistical analysis after normalisation samples using ANOVA+TUKEY ----
# 10.1 ANOVA to analyse whether there are differences between the groups
#     (Details see 4.) 

anova_pst_perc <- aov(prc_pos_avg ~ Pst_strain_n, data = df)

{summary(anova_pst_perc)
  summary.lm(anova_pst_perc)}

# 10.2 TUKEY test 
TUKEY_p <- TukeyHSD(anova_pst_perc)

generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  # put the labels in the same order as in the boxplot :
  Tukey.labels$Pst_strain = rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$Pst_strain) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset

LABELS_p <- generate_label_df(TUKEY_p , "Pst_strain_n")
names(LABELS_p) <- c("Letters_p","Pst_strain_n")
yvalue_p <- aggregate(Pst_CFU_n ~ Pst_strain_n, data = df, mean) # obtain letter position for y axis using means
LABELS_p$Pst_strain_n <- (unlist(lapply(strsplit(LABELS_p$Pst_strain_n,'[_]'),function(x){x[[1]]})))
final_p <- merge(yvalue_p,LABELS_p) 
final_p$Pst_CFU_n <- final_p$Pst_CFU_n + 1 # add +1 so that the label is not on the plot (e.g. outliers)
stat_labels_p <- as.vector(LABELS_n$Letters) # extract letter column as a vector to use for plot

# 10.3 safe TUKEY test results
sink(paste(pre, sep = "","CFU_count_DC3000_library_TUKEY_normalised_scaled.txt"))
print(LABELS_p)
sink()

# 11. EXPORT data for summary heat map ----

# 11.1 list all unique strains
strains <- unlist(unique(df$Pst_strain_n))

# 11.2 subset normalized simple_df by single vectors and calculate mean/sample

#strain <- "D36+HopA1"
#rm(ls=strain)
median_values <- lapply(strains, function(x){
data <- df[df$Pst_strain_n == x,]
medians <- median(data$prc_pos_avg)
return(medians)
})

# 11.3 generate small df to store information
strains <- as.character(strains)
df_median <- cbind("Pst_strain" = strains, "median" = median_values)

# 11.4 export if needed
write.csv(df_median, "./MM20230606_Median_CFU_counts_scaled.csv")

# 11.5 redo for scaled values



# a Dunnett's test after scaling
library("DescTools")
dt_no_anova <- DunnettTest(prc_pos_avg ~ Pst_strain_n, data = df, control = "D36", conf.level = 0.95)

# Plot the % data and see what changes

# 11 Box Plots comparing Pst strains, 2 dpi, no stats, normalised+scaled, save as svg (default) ----

# add objects to overlay horizontal lines on plot
D36E_df <- simple_df_norm[(simple_df_norm$Pst_strain_n == "D36"),]
D36E_median <- mean(D36E_df$Pst_CFU_n)
sd_plus <- D36E_median + sd(D36E_df$Pst_CFU_n)
sd_minus <- D36E_median - sd(D36E_df$Pst_CFU_n)

cbp1 <- c("#999999", "#44AA99", "#E69F00", "#661100", "#56B4E9", "#000000", "#117733", "#332288", "#F0E442", "#0072B2", "#D55E00", "#999933", "#CC79A7", "#E69F00")
cbp2 <- rep(x =  "#E69F00", times = 8)

# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# create pdf
#pdf(file = "C:/Users/Mende012/Documents/Bioinformatics/CFU_counts/MM20230107_CFU_count_norm_D36E.pdf",
#     width = 8,
#     height = 5)

# basic boxplot
g4 <-  ggplot(
  
  #input data
  df, aes(x=as.factor(Pst_strain_n), y=prc_pos_avg)) +
  
  # generate basic boxplot
  geom_boxplot(fill="white",                 # box colour
               outlier.colour = "white",     # Outliers color, 
               alpha=0) +                    # Box color transparency
  
  # overlay with jitter
  geom_jitter(shape=16, position=position_jitter(0.1),
              aes(colour = factor(Date))) +
  scale_colour_manual(values = cbp2) +
  
  # define the theme of the boxplot
  theme_bw() +  # make the bg white
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove background, frame
        axis.line = element_line(colour = "black")) +
  
  # label the axises 
  xlab("Pseudomonas strain") +                
  ylab("log10 fold change (CFU/cm2)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # turn the tick marks on the x axis 45 degree 
  

  # define axis limits if needed
  expand_limits( y = c(0, 5)) +
  
  # add statistical information: error bars, statistical results
  stat_boxplot(geom = "errorbar", # Error bars
                width = 0.2)  #+
  
  # ANOVA/TUKEY test labels
  #geom_text(data = final_n, aes(y=Pst_CFU_n, label = Letters_n), 
                              #vjust= -10 ,hjust= 0.5) 
  
  # significance for single comparisons
  #geom_signif(comparisons = list(c("D36", "D36+AvrE1"),
   #                              c("D36", "D36+HopM1"),
    #                             c("D36", "D36+HopU1"),
     #                            c("D36", "DC3000")),
      #        map_signif_level = TRUE,
       #       y_position = c(3, 3.5, 4, 4.5))+
  
              
  #  ggsignif::geom_signif(
  #  comparisons = list(c("D36", "D36+AvrE1"),
  #                     c("D36", "D36+HopM1"),
  #                     c("D36", "D36+HopU1"),
  #                     c("D36", "DC3000")),
  #  annotations = c("***", "***", "**", " ***"),
  #  y_position = c(3, 3.5, 4, 4.5) 
  
  #)
    
 # add ribbon to indicate st error of D36E
 # geom_ribbon(aes(ymin = sd_minus, ymax = sd_plus), alpha = 1) +  

# close pdf
ggsave(filename = paste(pre, sep = "","_CFU_count_normalised_scaled_nostats_mono.svg") ,width = 7, height = 4)
g4
dev.off()

ggplotly(g4)


# outlier detection and subsequent analysis
```



################################################################################
######################     DETECT and REMOVE OUTLIERS    #######################
################################################################################


```{r}
# calculate Q1 and Q3 for each treatment 
list_quantiles <- tapply(simple_df_norm$Pst_CFU_n, simple_df_norm$Pst_strain_n, quantile)
nr_samples <- length(list_quantiles)
 
Q1s <- sapply(1:nr_samples, function(i) list_quantiles[[i]][2])
Q3s <- sapply(1:nr_samples, function(i) list_quantiles[[i]][4])

# calculate IQR IQR = Q3 - Q1 for each treatement
IQRs <- tapply(simple_df_norm$Pst_CFU_n, simple_df_norm$Pst_strain_n, IQR)
  
# Calculate upper and lower limits
Lowers <- Q1s - 1.5*IQRs
Uppers <- Q3s + 1.5*IQRs

# generate df for upper and lower limits and merge them
df_lowers <- as.data.frame(Lowers)
df_uppers <- as.data.frame(Uppers)

df_lowers <- tibble::rownames_to_column(df_lowers, "Treatment")
df_uppers <- tibble::rownames_to_column(df_uppers, "Treatment")

IQR_Pst <- merge(df_lowers, df_uppers)
IQR_Pst$Treatment <- as.factor(IQR_Pst$Treatment) # convert Treatment into fctr so it matches the original 

# add IQUR upper/lower column to summary dataframe with all CFU values

df_outlier <- merge(simple_df_norm, IQR_Pst, by.x = "Pst_strain_n", by.y = "Treatment") #merge on bacteria

# new column for detecting outliers to
df_outlier$out <- ifelse(df_outlier$Pst_CFU_n < df_outlier$Lowers | df_outlier$Pst_CFU_n > df_outlier$Uppers, TRUE, FALSE )

# drop data outside the lower or upper range outside the IQR for each treatment 
df_no_outlier <- df_outlier[(df_outlier$out == FALSE),]

# plut the results after removing outliers 

# basic boxplot
g6 <-  ggplot(
  
  #input data
  df_no_outlier, aes(x=as.factor(Pst_strain_n), y=Pst_CFU_n)) +
  
  # generate basic boxplot
  geom_boxplot(fill="white",                 # box colour
               outlier.colour = "white",     # Outliers color, 
               alpha=0) +                    # Box color transparency
  
  # overlay with jitter
  geom_jitter(shape=16, position=position_jitter(0.1),
              aes(colour = factor(Date))) +
  scale_colour_manual(values = cbp1) +
  
  # define the theme of the boxplot
  theme_bw() +  # make the bg white
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove background, frame
        axis.line = element_line(colour = "black")) +
  
  # label the axises 
  xlab("Pseudomonas strain") +                
  ylab("log10 fold change (CFU/cm2)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # turn the tick marks on the x axis 45 degree 
  

  # define axis limits if needed
  expand_limits( y = c(0, 5)) +
  
  # add statistical information: error bars, statistical results
  stat_boxplot(geom = "errorbar", # Error bars
                width = 0.2)  +
  
  # ANOVA/TUKEY test labels
  #geom_text(data = final_n, aes(y=Pst_CFU_n, label = Letters_n), 
                              #vjust= -10 ,hjust= 0.5) 
  
  # significance for single comparisons
  #geom_signif(comparisons = list(c("D36", "D36+AvrE1"),
   #                              c("D36", "D36+HopM1"),
    #                             c("D36", "D36+HopU1"),
     #                            c("D36", "DC3000")),
      #        map_signif_level = TRUE,
       #       y_position = c(3, 3.5, 4, 4.5))+
  
              
    ggsignif::geom_signif(
    comparisons = list(c("D36", "D36+AvrE1"),
                       c("D36", "D36+HopM1"),
                       c("D36", "D36+HopU1"),
                       c("D36", "DC3000")),
    annotations = c("***", "***", "**", " ***"),
    y_position = c(3, 3.5, 4, 4.5) 
  
  )



# LAST: Safe stat results for better comparison in one df ----


stat_comp <- data.frame(  "Pst_strain" = LABELS$Pst_strain,
                          "not_normalised" = LABELS$Letters,
                          "normalised_neg_control" = LABELS_n$Letters_n,
                          "normalised_scaled" = LABELS_p$Letters)


sink(paste(pre, sep = "","CFU_count_DC3000_library_summary_stats.txt"))
print(stat_comp)
sink()




