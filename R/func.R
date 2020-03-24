
#################################################
# Functions used in Harrison et al. meta-analysis 
# 12 March 2020
#################################################

# This function will fit multi-level meta-analytic models (intercept only) to estimate the overall effect on the mean and CV accounting for study and species random effects and estimating a residual error (within study/species) variance. 

#TO DO: Incorporate the phylogeny when the phylogenies are corrected and available
fit_int_MLMAmodels <- function(data, phylo_vcv){

          lnCVR <- metafor::rma.mv(CVR_yi ~ 1, V = CVR_vi, random = list(~1|study_ID, ~1|spp, ~1|obs), R = list(spp=phylo_vcv), data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ 1, V = SMD_vi, random = list(~1|study_ID, ~1|spp, ~1|obs), R = list(spp=phylo_vcv), data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait <- function(data){
          lnCVR <- metafor::rma.mv(CVR_yi ~ personality_trait, V = CVR_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ personality_trait, V = SMD_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait_SSD_index <- function(data){
          lnCVR <- metafor::rma.mv(CVR_yi ~ personality_trait*SSD_index, V = CVR_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ personality_trait*SSD_index, V = SMD_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait_mating_system  <- function(data){
          lnCVR <- metafor::rma.mv(CVR_yi ~ personality_trait*mating_system, V = CVR_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ personality_trait*mating_system, V = SMD_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_parental_care  <- function(data){
          lnCVR <- metafor::rma.mv(CVR_yi ~ parental_care, V = CVR_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ parental_care, V = SMD_vi, random = list(~1|study_ID, ~1|species_name, ~1|obs), data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

# This function takes the dataset, splits the dataset up by broad taxonomic group and applies a specified meta-analytic or meta-regression model to each taxa, saving the model results into a list that can be used down stream
meta_model_fits <- function(data, phylo_vcv, type = c("int", "pers", "pers_SSD", "pers_mate", "parent_care")){
         taxa_list <- split(data, data$taxo_group)

      type <-  match.arg(type)

    if(type == "int"){
        model_fits <- mapply(function(x, y) fit_int_MLMAmodels(x, y), x = taxa_list, y = phylo_vcv)
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

tree_checks <- function(data, tree){
        # How many unqiue species exist in data and tree
        Numbers <- matrix(,nrow = 2, ncol = 2)
        Numbers[1,2] <- length(unique(data$spp_name_phylo)) #108
        Numbers[2,2] <- length(tree$tip.label) #109

        Numbers[,1] <- c("Species in data:", "Species in tree:")

        # Missing species or species not spelt correct
        #species_list <- tree$tip.label[!tree$tip.label %in% unique(data$spp)]        
        species_list1= setdiff(sort(tree$tip.label), sort(unique(data$spp_name_phylo)))
        species_list2= setdiff(sort(unique(data$spp_name_phylo)), sort(tree$tip.label) )

        return(list(SpeciesNumbers = Numbers, SpeciesList_NotFound_InData_But_NotTree=species_list1, SpeciesList_NotFound_InTree_But_NotData=species_list2))
      }


    read_birds <- function(x){
          tr <- read.nexus(x)
          ave.tree <- midpoint(ls.consensus(tr))
          return(ave.tree)
    }
