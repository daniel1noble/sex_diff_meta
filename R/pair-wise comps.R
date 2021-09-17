#### Running pairwise comparisons on meta-analysis models ####

# Here we are comparing the mean lnCVR estimates for bird and mammal personalities to compare the effect of heterogamety.
# We want to compare boldness in birds to boldness in mammals, etc. to look for an effect of heterogamety on personality variance.
# Since these were two seperate models we have to do a bit of adjusting in order to run pairwise comparisons. 

# Option 1 
# hand calculate contrasts 
  # function to calculate contrasts, just feed in the relevant data from model output
  contrast_means <- function(m1, m2, L95m1, U95m1, L95m2, U95m2){
    # Calculate contrast
    diff = (m2 - m1)
  
    # Calculate SE from CI
    se_m1 = ((U95m1 - L95m1)) / (2*1.96) 
    se_m2 = ((U95m2 - L95m2)) / (2*1.96)
  
    # Calculate SE for contrast
    se_diff = sqrt(se_m1^2 + se_m2^2)
  
    return(data.frame(contrast = diff, se = se_diff))
}

contrasts <- contrast_means(m1 =-0.3095, m2 =0.0782, L95m1 =-0.7558, U95m1 =0.1368, L95m2=-0.5244, U95m2=0.6808)

  # If you want inferential stats from the
  p = with(contrasts, dt(contrast/se, df = 1154))

# Option 2 - formal pairwise comparison model
# we use the glht function in the multcomp package to run Tukey tests

library(multcomp)
library(Matrix)

# Our phylogenies were seperate for each taxo group so need to combine birds and mammals into 1 Matrix
  # Create block diag phylogeny
  phylogeny <- Matrix::bdiag(phylo_vcv_bird, phylo_vcv_mammal) # use this as the phylo vcv

  # needs to have colnames for use in random effects model
  dimnames(phylogeny) <- Map(c, dimnames(phylo_vcv_bird), dimnames(phylo_vcv_mammal))

# Creating the model 
# There should be no intercept because we want to compare the means directly, with pers trait and taxo group as mods (including their interaction)

  # only include bird and mammal data
  pers_new_contrast <- as.data.frame(pers_new %>%
                                     filter(taxo_group =="mammal" | taxo_group == "bird"))

  # lnCVR model
  contrast_birdmammal_lncvr <- rma.mv(CVR_yi ~ personality_trait*taxo_group -1, V = CVR_vi, 
                                    random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), 
                                    R = list(spp_name_phylo=phylogeny), control=list(optimizer="optim"), 
                                    test = "t", data = pers_new_contrast)

  # Example of setting up contrasts of interest. 
  summary(glht(contrast_birdmammal_lncvr, linfct=rbind(b1=c(1,0,0,0,0,0,0,0,0,0), b2=c(0,0,0,0,0,1,0,0,0,0))), test=Chisqtest())

  # our comparison
  summary(glht(contrast_birdmammal_lncvr, linfct = cbind(contrMat(rep(1,10), type = "Tukey"))), test=adjusted("none"))
  
  contrast_heterogamety <- glht(contrast_birdmammal_lncvr, linfct=mcp(TimeGeno="Tukey"))

# The output from our hand-calculations should roughly line up with the more formal comparison but they won't match exactly 
# because it does not take into account random effects. It is just a test to make sure we are comparing the mean differences
