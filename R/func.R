
#################################################
# Functions used in Harrison et al. meta-analysis 
# 12 March 2020
#################################################

# This function will fit multi-level meta-analytic models (intercept only) to estimate the overall effect on the mean and CV accounting for study and species random effects and estimating a residual error (within study/species) variance. 

#TO DO: Incorporate the phylogeny when the phylogenies are corrected and available
fit_int_MLMAmodels <- function(data){
          lnCVR <- metafor::rma.mv(CVR_yi ~ 1, V = CVR_vi, random = list(~1|study_ID, ~1|species_name, ~1|err), data = data)

            SMD <- metafor::rma.mv(SMD_yi ~ 1, V = SMD_vi, random = list(~1|study_ID, ~1|species_name, ~1|err), data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}


# This function takes the dataset, splits the dataset up by broad taxonomic group and applies a specified meta-analytic or meta-regression model to each taxa, saving the model results into a list that can be used down stream
meta_model_fits <- function(data){
         taxa_list <- split(pers_new, pers_new$taxo_group.x)

        model_fits <- lapply(taxa_list, function(x) fit_int_MLMAmodels(x))
} 


