
#################################################
# Functions used in Harrison et al. meta-analysis 
# 12 March 2020
#################################################

# Theses functions will fit multi-level meta-analytic models (intercept only) and meta-regression models to estimate the overall effect on the mean and CV accounting for study and species (phylogeny) random effects and estimating a residual error (within study/species) variance. 

fit_int_MLMAmodels <- function(data, phylo_vcv){

          lnCVR <- metafor::rma.mv(CVR_yi ~ 1, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ 1, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait <- function(data, phylo_vcv){
          lnCVR <- metafor::rma.mv(CVR_yi ~ personality_trait-1, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ personality_trait-1, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait_SSD_index <- function(data, phylo_vcv){
          lnCVR <- metafor::rma.mv(CVR_yi ~ -1+ personality_trait + SSD_index, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ -1 + personality_trait + SSD_index, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait_mating_system  <- function(data, phylo_vcv){
          lnCVR <- metafor::rma.mv(CVR_yi ~ -1 + personality_trait + mating_system, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ -1 + personality_trait + mating_system, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_parental_care  <- function(data, phylo_vcv){
          lnCVR <- metafor::rma.mv(CVR_yi ~ parental_care - 1, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ parental_care - 1, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

sensitivity_models <- function(){
  
}

# This function takes the dataset, splits the dataset up by broad taxonomic group and applies a specified meta-analytic or meta-regression model to each taxa (seee functions for models above), saving the model results into a list that can be used down stream
meta_model_fits <- function(data, phylo_vcv, type = c("int", "pers", "pers_SSD", "pers_mate", "parent_care")){
         taxa_list <- split(data, data$taxo_group)

      type <-  match.arg(type)

    if(type == "int"){
        model_fits <- mapply(function(x, y) fit_int_MLMAmodels(x, y), x = taxa_list, y = phylo_vcv)
    }

    if(type == "pers"){
        model_fits <- mapply(function(x, y) fit_MLMA_reg_models_personality_trait(x, y), x = taxa_list, y = phylo_vcv)
    }

  	if(type == "pers_SSD"){
        model_fits <- mapply(function(x,y) fit_MLMA_reg_models_personality_trait_SSD_index(x, y), x = taxa_list, y = phylo_vcv)
    }

	if(type == "pers_mate"){
        model_fits <- mapply(function(x,y) fit_MLMA_reg_models_personality_trait_mating_system(x, y), x = taxa_list, y = phylo_vcv)
    }    

    if(type == "parent_care"){
        model_fits <- mapply(function(x, y) fit_MLMA_reg_models_parental_care(x, y), x = taxa_list, y = phylo_vcv)
    }  

   return(model_fits)
} 

# This function is very useful for both checking that taxa names in data and tree match, printing out where discrepancies lie and then allowing the trees to be pruned to the taxa within the dataset assuming all other taxa match. 

tree_checks <- function(data, tree, dataCol, type = c("checks", "prune")){
        type = match.arg(type)

        # How many unique species exist in data and tree
        Numbers <- matrix(nrow = 2, ncol = 1)
        Numbers[1,1] <- length(unique(data[,dataCol])) 
        Numbers[2,1] <- length(tree$tip.label) 

        rownames(Numbers)<- c("Species in data:", "Species in tree:")

        # Missing species or species not spelt correct      
        species_list1= setdiff(sort(tree$tip.label), sort(unique(data[,dataCol])))
        species_list2= setdiff(sort(unique(data[,dataCol])), sort(tree$tip.label) )

        if(type == "checks"){
          return(list(SpeciesNumbers = data.frame(Numbers), 
                    Species_InTree_But_NotData=species_list1, 
                    Species_InData_But_NotTree=species_list2))
          }
        
        if(type == "prune"){
            if(length(species_list2) >=1) stop("Sorry, you can only prune a tree when you have no taxa existing in the data that are not in the tree")

              return(ape::drop.tip(tree, species_list1))
              }
      }



# Function that will read in nexus file of 100 bird trees and take a consensus / average bird tree from these
    read_birds <- function(x){
          tr <- read.nexus(x)
          ave.tree <- midpoint(ls.consensus(tr))
          return(ave.tree)
    }
