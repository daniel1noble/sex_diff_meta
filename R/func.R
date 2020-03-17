
#################################################
# Functions used in Harrison et al. meta-analysis 
# 12 March 2020
#################################################

# This function will fit multi-level meta-analytic models (intercept only) to estimate the overall effect on the mean and CV accounting for study and species random effects and estimating a residual error (within study/species) variance. 

#TO DO: Incorporate the phylogeny when the phylogenies are corrected and available
fit_int_MLMAmodels <- function(data){
          lnCVR <- metafor::rma.mv(CVR_yi ~ 1, V = CVR_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data)

            SMD <- metafor::rma.mv(SMD_yi ~ 1, V = SMD_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait <- function(data){
          lnCVR <- metafor::rma.mv(CVR_yi ~ personality_trait, V = CVR_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data)

            SMD <- metafor::rma.mv(SMD_yi ~ personality_trait, V = SMD_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait_SSD_index <- function(data){
          lnCVR <- metafor::rma.mv(CVR_yi ~ personality_trait*SSD_index, V = CVR_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data)

            SMD <- metafor::rma.mv(SMD_yi ~ personality_trait*SSD_index, V = SMD_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait_mating_system  <- function(data){
          lnCVR <- metafor::rma.mv(CVR_yi ~ personality_trait*mating_system, V = CVR_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data)

            SMD <- metafor::rma.mv(SMD_yi ~ personality_trait*mating_system, V = SMD_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_parental_care  <- function(data){
          lnCVR <- metafor::rma.mv(CVR_yi ~ parental_care, V = CVR_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data)

            SMD <- metafor::rma.mv(SMD_yi ~ parental_care, V = SMD_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

# This function takes the dataset, splits the dataset up by broad taxonomic group and applies a specified meta-analytic or meta-regression model to each taxa, saving the model results into a list that can be used down stream
meta_model_fits <- function(data, type = c("int", "pers", "pers_SSD", "pers_mate", "parent_care")){
         taxa_list <- split(pers_new, pers_new$taxo_group.x)

      tyep <-  match.arg(type)

    if(type == "int"){
        model_fits <- lapply(taxa_list, function(x) fit_int_MLMAmodels(x))
    }

    if(type == "pers"){
        model_fits <- lapply(taxa_list, function(x) fit_MLMA_reg_models_personality_trait(x))
    }

  	if(type == "pers_SSD"){
        model_fits <- lapply(taxa_list, function(x) fit_MLMA_reg_models_personality_trait_SSD_index(x))
    }

	if(type == "pers_mate"){
        model_fits <- lapply(taxa_list, function(x) fit_MLMA_reg_models_personality_trait_mating_system(x))
    }    

    if(type == "parent_care"){
        model_fits <- lapply(taxa_list, function(x) fit_MLMA_reg_models_parental_care(x))
    }  

   return(model_fits)
} 


