#################################################
# Setting up dataset used in Harrison et al. personality meta-analysis 
# 31 August 2020
#################################################

# These datasets and the following code are used to create the final dataset 'pers_new' used for further analysis. 

#################################################

# Clear work space
rm(list=ls())

# Load packages
install.packages("pacman")
pacman::p_load(knitr, dplyr, kableExtra, tidyverse, metafor, ggplot2)

# Source our own functions
source("./R/func.R")

## 1. setting up dataset

  # Import original datasets 
  pers <- read.csv("./data/pers_data.csv", stringsAsFactors = FALSE)
  bodysize <- read.csv("./data/bodysize_SSD.csv", stringsAsFactors = FALSE)

  # Merge the two by spp_names columns
  pers <- merge(x = pers,
              y = bodysize[,c("species_name", "SSD_index", "mating_system")],
              by="species_name", all.x=TRUE,  no.dups = TRUE)

  # Select the relevant columns to make things easier
  pers_new <- pers %>% 
  select(study_ID, year, species_name, SSD_index, taxo_group, data_type, personality_trait, male_n, male_mean_conv, 
         male_sd_conv, female_n, female_mean_conv, female_sd_conv, depend, directionality, spp_name_phylo, mating_system, 
         age, population, study_environment, study_type, measurement_type)

## 2. Look at dataset to make sure everything is looking good

  # Check species numbers
  pers_new %>%
  group_by(taxo_group) %>%
  summarise(species = length(unique(species_name)))

  # Add in observation level random effect (metafor doesn't do this, need to do it manually)
  pers_new <- pers_new %>% 
  group_by(taxo_group) %>% 
  mutate(obs = 1:length(study_ID))

## ISSUE: 
# Our dataset has 1) latency data which has not been corrected for normality (i.e. log-transformed), and 
# 2) proportional data which constrains distributions of variances at high and low values. 
# As such, we need to correct these means and SDs in the dataset before we calculate SMD and lnCVR effect sizes. 
# Scores will be considered as 'normal' data. These corrections should get rid of any outliers.

  # converting proportions to logit scale 
    # filter dataset to show only proportional data
    pers_proportion <- as.data.frame(pers_new %>%
                                   filter(measurement_type == "proportion"))

    # both sexes combined in the one function now - save list to object so I can save as dataframe
    proportion_all <- convert_propor(pers_proportion$male_mean, pers_proportion$m_SE_to_SD, 
                                 pers_proportion$female_mean, pers_proportion$f_SE_to_SD)

  # converting latency data to log-latency
    # filter dataset to show only latency data
    pers_latency <- as.data.frame(pers_new %>%
                                filter(measurement_type == "latency"))

    # both sexes combined in the one function now - save list to object so we can save as dataframe
    latency_all <- convert_latency(pers_latency$male_mean, pers_latency$m_SE_to_SD, pers_latency$female_mean, 
                               pers_latency$f_SE_to_SD, method = "analyt") 

    # a bunch of NAs for some reason, but running the calculation with the means and SD works 
    # NOTE did the NAs manually by plugging mean and SD into calculation provided in func.R 

  # make table into dataframe, then add in manually
  write.csv(proportion_all, file = "pers_proportion.csv", row.names = FALSE)  
  write.csv(latency_all, file = "pers_latency.csv", row.names = FALSE)  
  
## ISSUE: 
# We have 3 papers in our eligible studies dataset (P172, P210 and P231) that might have issues with data duplication / unreliability
# LMH checked the retraction database retractiondatabase.org frequently to check if these were retracted or flagged.
# As of 31 August 2020 none of the 3 papers had been retracted, so we decided to keep these studies in our final dataset. 

## 3. calculating effect sizes (SMD & lnCVR)  
  
  # SMD (Hedge's g)
  pers_new <- escalc(measure = "SMD", 
                     n1i = male_n, n2i = female_n,
                     m1i = male_mean_conv, m2i = female_mean_conv,
                     sd1i = male_sd_conv, sd2i = female_sd_conv, data = pers_new, var.names=c("SMD_yi","SMD_vi"), append = TRUE)
  
  # lnCVR
  pers_new <- escalc(measure = "CVR",
                     n2i = female_n, n1i = male_n,
                     m2i = female_mean_conv, m1i = male_mean_conv,
                     sd2i = female_sd_conv, sd1i = male_sd_conv, data = pers_new, var.names=c("CVR_yi","CVR_vi"))

  # we have some NAs where one or both sexes have a value of 0 for either mean or SD. Will be easiest to just remove these.
    # Exclude NAs
    pers_new <- pers_new %>%
    filter(!is.na(CVR_yi), !is.na(SMD_yi))
    dim(pers_new)    # check they've been removed with no issues
    
## 4. Calculating the mean-variance relationship
# This is to make sure that our using lnCVR as a measure of variability is valid 
    # females and males seperately because they are in different columns
    # use ggplot to make a scatterplot of females
    fem <- ggplot(pers_new, aes(x = female_mean_conv, y = female_sd_conv)) + geom_point()
    
    # on log scale
    fem + scale_x_continuous(trans = 'log10')
      + scale_y_continuous(trans = 'log10')
    
    # mean and SD on log scale to calculate correlation
    logfemale_mean <- log(pers_new$female_mean_conv)
    logfemale_SD <- log(pers_new$female_sd_conv)
    
    # correlation between mean and SD
    cor(logfemale_mean, logfemale_SD) #0.91
  
  # Males
    # use ggplot to make a scatterplot of females
    male <- ggplot(pers_new, aes(x = male_mean_conv, y = male_sd_conv)) + geom_point()
    
    # on log scale
    male + scale_x_continuous(trans = 'log10')
    + scale_y_continuous(trans = 'log10')
    
    # mean and SD on log scale to calculate correlation
    logmale_mean <- log(pers_new$male_mean_conv)
    logmale_SD <- log(pers_new$male_sd_conv)
    
    # correlation between mean and SD
    cor(logmale_mean, logmale_SD) #0.90

## 5. Checking for outliers
  # lnCVR
    hist(pers_new$CVR_yi)
    funnel(x = pers_new$CVR_yi, vi = pers_new$CVR_vi, yaxis="seinv", xlim = c(-10, 10))
  
    #SMD
    hist(pers_new$SMD_yi)
    funnel(x = pers_new$SMD_yi, vi = pers_new$SMD_vi, yaxis="seinv")
    text(as.character(pers_new$study_ID), x = pers_new$SMD_yi, y = 1/sqrt(pers_new$SMD_vi))
    
    
    pers_new %>% 
      filter(study_ID == "P231" & SMD_yi < -20) # P231 has 2 datapoints (proportions) with huge differences: -68 and -38 remove?
    
    # filter out those two large effect sizes
    pers_new <- pers_new %>%
      filter(!SMD_yi < -20)
    dim(pers_new) 
    
    # check funnel once more 
    pers_new %>% 
      filter(study_ID == "P206" & SMD_yi < -15) # looks like a genuine effect size, I'm going to leave it in
  
    # checking SMD outliers - inverse SE > 14
    funnel(x = pers_new$SMD_yi, vi = pers_new$SMD_vi, yaxis="seinv")
    text(as.character(pers_new$obs), x = pers_new$SMD_yi, y = 1/sqrt(pers_new$SMD_vi), offset = 0.8)
    
    # checking lnCVR outliers
    funnel(x = pers_new$CVR_yi, vi = pers_new$CVR_vi, yaxis="seinv", xlim = c(-10, 10))
    text(as.character(pers_new$obs), x = pers_new$CVR_yi, y = 1/sqrt(pers_new$CVR_vi), offset = 0.8) 
      
 # I noticed some measures are more physiological/not personality than personality, so I think it's wise to remove these before we run the models - just to be safe!
    # P029 - obs 22, 23, 32 
    # P084 - obs 59, 62, 63, 65, 68, 70, 71, 72, 74
    # P060 - obs 216, 217
    # P211 - obs 230, 245
    # P117 - obs 397, 393, 402, 414
    # P197 - obs 541, 544, 546, 547
    # P069 - obs 669, 672, 673, 682, 683, 684, 686, 694 
    
    # remove these by study 
    pers_new <- pers_new %>% filter(!study_ID == "P029" | !obs %in% c("22", "23", "32"))
    
    pers_new <- pers_new %>% filter(!study_ID == "P084" | !obs %in% c("59", "62", "63", "65", "68", "70", "71", "72", "74"))
    
    pers_new <-  pers_new %>% filter(!study_ID == "P060" | !obs %in% c("216", "217"))
    
    pers_new <- pers_new %>% filter(!study_ID == "P211" | !obs %in% c("230", "245"))
    
    pers_new <- pers_new %>% filter(!study_ID == "P117" | !obs %in% c("397", "393", "402", "414"))
    
    pers_new <- pers_new %>% filter(!study_ID == "P197" | !obs %in% c("541", "544", "546", "547"))
    
    pers_new <- pers_new %>% filter(!study_ID == "P197" | !obs %in% c("669", "672", "673", "682", "683", "684", "686", "694")) 
    
    dim(pers_new) 
        
# see 'pers_new.csv' for final dataset after all these checks

