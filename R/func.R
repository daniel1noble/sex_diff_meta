
#################################################
# Functions used in Harrison et al. meta-analysis 
# 12 March 2020
#################################################

# These functions will fit multi-level meta-analytic models (intercept only) and meta-regression models to estimate the overall effect on the mean and CV accounting for study and species (phylogeny) random effects and estimating a residual error (within study/species) variance. 

fit_int_MLMAmodels <- function(data, phylo_vcv){

          lnCVR <- metafor::rma.mv(CVR_yi ~ 1, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ 1, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait <- function(data, phylo_vcv){
          lnCVR <- metafor::rma.mv(CVR_yi ~ personality_trait-1, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ personality_trait-1, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_pubbias <- function(data, phylo_vcv){
  lnCVR <- metafor::rma.mv(CVR_yi ~ -1 + personality_trait + sqrt_inv_n_tilda, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data)
  
  SMD <- metafor::rma.mv(SMD_yi_flip ~ -1 + personality_trait + sqrt_inv_n_tilda, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data) 
  
  return(list(SMD = SMD, 
              lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait_SSD_index <- function(data, phylo_vcv){
          lnCVR <- metafor::rma.mv(CVR_yi ~ personality_trait*SSD_index, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ personality_trait*SSD_index, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait_SSD_simple <- function(data, phylo_vcv){
  lnCVR <- metafor::rma.mv(CVR_yi ~ SSD_index, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data)
  
  SMD <- metafor::rma.mv(SMD_yi_flip ~ SSD_index, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data) 
  
  return(list(SMD = SMD, 
              lnCVR = lnCVR))
}

fit_MLMA_reg_models_mating_system  <- function(data, phylo_vcv){
          lnCVR <- metafor::rma.mv(CVR_yi ~ mating_system, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ mating_system, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_parental_care  <- function(data, phylo_vcv){
          lnCVR <- metafor::rma.mv(CVR_yi ~ parental_care - 1, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data)

            SMD <- metafor::rma.mv(SMD_yi_flip ~ parental_care - 1, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), test = "t", data = data) 

          return(list(SMD = SMD, 
          			lnCVR = lnCVR))
}

fit_MLMA_reg_models_age  <- function(data, phylo_vcv){
  lnCVR <- metafor::rma.mv(CVR_yi ~ age, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data)
  
  SMD <- metafor::rma.mv(SMD_yi_flip ~ age, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data) 
  
  return(list(SMD = SMD, 
              lnCVR = lnCVR))
}

fit_MLMA_reg_models_population <- function(data, phylo_vcv){
  lnCVR <- metafor::rma.mv(CVR_yi ~ population, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data)
  
  SMD <- metafor::rma.mv(SMD_yi_flip ~ population, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data) 
  
  return(list(SMD = SMD, 
              lnCVR = lnCVR))
}

fit_MLMA_reg_models_environment <- function(data, phylo_vcv){
  lnCVR <- metafor::rma.mv(CVR_yi ~ study_environment, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data)
  
  SMD <- metafor::rma.mv(SMD_yi_flip ~ study_environment, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data) 
  
  return(list(SMD = SMD, 
              lnCVR = lnCVR))
}

fit_MLMA_reg_models_studytype <- function(data, phylo_vcv){
  lnCVR <- metafor::rma.mv(CVR_yi ~ study_type, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data)
  
  SMD <- metafor::rma.mv(SMD_yi_flip ~ study_type, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data) 
  
  return(list(SMD = SMD, 
              lnCVR = lnCVR))
}

fit_int_MLMAmodel_D <- function(data, phylo_vcv, D){
    lnCVR <- mapply(function(x,y,z) metafor::rma.mv(CVR_yi ~ 1, V = CVR_vi, 
                                                  random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), 
                                                  R = list(spp_name_phylo=y, obs = z), 
                                                  data = x), x = split_taxa, y = phylo_vcv, z = D)
  
  SMD <- mapply(function(x,y,z) metafor::rma.mv(SMD_yi_flip ~ 1, V = SMD_vi, 
                                                random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), 
                                                R = list(spp_name_phylo=y, obs = z), 
                                                data = x), x = split_taxa, y = phylo_vcv, z = D) 
  
  return(list(SMD = SMD, 
              lnCVR = lnCVR))
}

fit_MLMA_reg_models_personality_trait_D <- function(data, phylo_vcv, D){
                lnCVR <- mapply(function(x,y,z) metafor::rma.mv(CVR_yi ~ personality_trait-1, V = CVR_vi, 
                                                  random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), 
                                                  R = list(spp_name_phylo=y, obs = z), 
                                                  data = x), x = split_taxa, y = phylo_vcv, z = D_matrices_0.8)
  
                SMD <- mapply(function(x,y,z) metafor::rma.mv(SMD_yi_flip ~ personality_trait-1, V = SMD_vi, 
                                                random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), 
                                                R = list(spp_name_phylo=y, obs = z), 
                                                data = x), x = split_taxa, y = phylo_vcv, z = D_matrices_0.8) 
  
        return(list(SMD = SMD, 
              lnCVR = lnCVR))
}

sensitivity_mod1 <- function(data, phylo_vcv){
  lnCVR <- metafor::rma.mv(CVR_yi ~ score, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data)
  
  SMD <- metafor::rma.mv(SMD_yi_flip ~ score, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data) 
  
  return(list(SMD = SMD, 
              lnCVR = lnCVR))
}

sensitivity_mod2 <- function(data, phylo_vcv){
  lnCVR <- metafor::rma.mv(CVR_yi ~ proportion, V = CVR_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data)
  
  SMD <- metafor::rma.mv(SMD_yi_flip ~ proportion, V = SMD_vi, random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), R = list(spp_name_phylo=phylo_vcv), control=list(optimizer="optim"), test = "t", data = data) 
  
  return(list(SMD = SMD, 
              lnCVR = lnCVR))
}

# This function takes the dataset, splits the dataset up by broad taxonomic group and applies a specified meta-analytic or meta-regression model to each taxa (seee functions for models above), saving the model results into a list that can be used down stream
meta_model_fits <- function(data, phylo_vcv, type = c("int", "pers", "pers_SSD", "pers_mate", "parent_care", "age", "pop", "environ", "study_type", "score", "proportion", "pubbias")){
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
        model_fits <- mapply(function(x,y) fit_MLMA_reg_models_mating_system(x, y), x = taxa_list, y = phylo_vcv)
    }    

    if(type == "parent_care"){
        model_fits <- mapply(function(x, y) fit_MLMA_reg_models_parental_care(x, y), x = taxa_list, y = phylo_vcv)
    }  

      if(type == "age"){
        model_fits <- mapply(function(x, y) fit_MLMA_reg_models_age(x, y), x = taxa_list, y = phylo_vcv)
      }    
     
      if(type == "pop"){
        model_fits <- mapply(function(x, y) fit_MLMA_reg_models_population(x, y), x = taxa_list, y = phylo_vcv)
      }  
      
      if(type == "environ"){
        model_fits <- mapply(function(x, y) fit_MLMA_reg_models_environment(x, y), x = taxa_list, y = phylo_vcv)
      }  
      
      if(type == "study_type"){
        model_fits <- mapply(function(x, y) fit_MLMA_reg_models_studytype(x, y), x = taxa_list, y = phylo_vcv)
      }  
      
      if(type == "score"){
        model_fits <- mapply(function(x, y) sensitivity_mod1(x, y), x = taxa_list, y = phylo_vcv)
      } 
      
      if(type == "proportion"){
        model_fits <- mapply(function(x, y) sensitivity_mod2(x, y), x = taxa_list, y = phylo_vcv)
      } 
      
      if(type == "pubbias"){
        model_fits <- mapply(function(x, y) fit_MLMA_reg_models_pubbias(x, y), x = taxa_list, y = phylo_vcv)
      } 
      
   return(model_fits)
} 

# D matrix intercept models
fit_int_MLMAmodel_D <- function(data, phylo_vcv, D){
  
  lnCVR <- mapply(function(x,y,z) metafor::rma.mv(CVR_yi ~ 1, V = CVR_vi, 
                                                  random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), 
                                                  R = list(spp_name_phylo=y, obs = z), 
                                                  data = x, test = "t"), x = split_taxa, y = phylo_vcv, z = D, SIMPLIFY = FALSE)
  
  SMD <- mapply(function(x,y,z) metafor::rma.mv(SMD_yi_flip ~ 1, V = SMD_vi, 
                                                random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), 
                                                R = list(spp_name_phylo=y, obs = z), 
                                                data = x, test = "t"), x = split_taxa, y = phylo_vcv, z = D, SIMPLIFY = FALSE) 
  
  return(list(SMD=SMD, 
              lnCVR = lnCVR))
}

# D matrix personality trait models
fit_int_MLMAmodel_D_pers <- function(data, phylo_vcv, D){
  
  lnCVR <- mapply(function(x,y,z) metafor::rma.mv(CVR_yi ~ personality_trait-1, V = CVR_vi, 
                                                  random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), 
                                                  R = list(spp_name_phylo=y, obs = z), 
                                                  data = x, test = "t"), x = split_taxa, y = phylo_vcv, z = D, SIMPLIFY = FALSE)
  
  SMD <- mapply(function(x,y,z) metafor::rma.mv(SMD_yi_flip ~ personality_trait-1, V = SMD_vi, 
                                                random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), 
                                                R = list(spp_name_phylo=y, obs = z), 
                                                data = x, test = "t"), x = split_taxa, y = phylo_vcv, z = D, SIMPLIFY = FALSE) 
  
  return(list(SMD=SMD, 
              lnCVR = lnCVR))
}

# D matrix personality trait x SSD models
fit_int_MLMAmodel_D_pers_ssd <- function(data, phylo_vcv, D){
  
  lnCVR <- mapply(function(x,y,z) metafor::rma.mv(CVR_yi ~ personality_trait*SSD_index, V = CVR_vi, 
                                                  random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), 
                                                  R = list(spp_name_phylo=y, obs = z), control=list(optimizer="optim"),
                                                  data = x, test = "t"), x = split_taxa, y = phylo_vcv, z = D, SIMPLIFY = FALSE)
  
  SMD <- mapply(function(x,y,z) metafor::rma.mv(SMD_yi_flip ~ personality_trait*SSD_index, V = SMD_vi, 
                                                random = list(~1|study_ID, ~1|spp_name_phylo, ~1|obs), 
                                                R = list(spp_name_phylo=y, obs = z), control=list(optimizer="optim"),
                                                data = x, test = "t"), x = split_taxa, y = phylo_vcv, z = D, SIMPLIFY = FALSE) 
  
  return(list(SMD=SMD, 
              lnCVR = lnCVR))
}

meta_model_fits_D <- function(data, phylo_vcv, D, type = c("int", "pers", "pers_SSD")){
  taxa_list <- split(data, data$taxo_group)
  
  type <-  match.arg(type)
  
  if(type == "int"){
    model_fits_D <- mapply(function(x, y, z) fit_int_MLMAmodel_D(x, y, z), x = taxa_list, y = phylo_vcv, z = D)
  }
  
  if(type == "pers"){
    model_fits_D <- mapply(function(x, y, z) fit_MLMA_reg_models_personality_trait_D(x, y, z), x = taxa_list, y = phylo_vcv, z = D)
  }
  
  if(type == "pers_SSD"){
    model_fits_D <- mapply(function(x, y, z) fit_MLMA_reg_models_personality_trait_SSD_index_D(x, y, z), x = taxa_list, y = phylo_vcv, z = D)
  }
  
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

# modified orchard function to get pred intervals from MLMR model with continuous moderator term:
    # Get prediction intervals for a single categorical variable when controlling for continuous variables in a model.
    get_pred_est_cont <- function (model, mod_cat, mod_cont, type = c("mean", "zero")) {
      type = match.arg(type)
      # Categorical moderator
      pos <- str_which(row.names(model$beta), {{mod_cat}}, negate = FALSE)
      name_cat <- row.names(model$beta)[pos]
      name_cat <- as.character(stringr::str_replace(name_cat, {{mod_cat}}, ""))
      name_clean <- firstup(name_cat)
      len <- length(name_clean)
      
      # Continuous moderators
      len_cont <- length(mod_cont)
      
      if(type == "mean"){
        if(len_cont == 1){  	
          means_cont <- mean(model$X[,mod_cont])
          names(means_cont) <- mod_cont
        }
        if(len_cont > 1){
          means_cont <- colMeans(model$X[,mod_cont])
        }
      }
      
      if(type == "zero"){
        if(len_cont == 1){  	
          means_cont <- 0
          names(means_cont) <- mod_cont
        }
        if(len_cont > 1){
          means_cont <- rep(0, len_cont)
          names(means_cont) <- mod_cont
        }
      }
      
      # Get prediction intervals
      newdata <- matrix(0, ncol = len+len_cont, nrow = len)
      diag(newdata[c(1:len),c(1:len)]) <- 1
      newdata[,c((len+1):(len+len_cont))] <- means_cont
      
      for(i in 1:len){
        newdata[i,c((len+1):(len+len_cont))] <-  means_cont
      }
      
      colnames(newdata) <- c(name_cat, names(means_cont))
      
      pred <- metafor::predict.rma(model, newmods = newdata)
      
      estimate <- pred$pred
      lowerCL <- pred$ci.lb 
      upperCL <- pred$ci.ub
      lowerPR <- pred$cr.lb
      upperPR <- pred$cr.ub 
      
      table <- tibble::tibble(name = name_clean, estimate = estimate, lowerCL = lowerCL, upperCL = upperCL, lowerPR = lowerPR, upperPR = upperPR)
      return(table)
    }
    
    #get_pred_est_cont(model, mod_cat, mod_cont)
    
    get_data_cont <- function(model, mod_cat){
      
      pos <- str_which(row.names(model$beta), {{mod_cat}}, negate = FALSE)
      names <- row.names(model$beta)[pos]
      name_cat <- as.character(stringr::str_replace(names, {{mod_cat}}, ""))
      
      X <- as.data.frame(model$X)[,names]
      
      moderator <- matrix(ncol = 1, nrow = dim(X)[1])
      
      for(i in 1:ncol(X)){
        moderator <- ifelse(X[,i] == 1, name_cat[i], moderator)
      }
      moderator <- firstup(moderator)
      yi <- model$yi
      vi <- model$vi
      type <- attr(model$yi, "measure")
      
      data <- data.frame(yi, vi, moderator, type)
      return(data)
      
    }
    
    mod_results_new <- function(model, mod_cat, mod_cont, type = c("mean", "zero")) { 
      
      type = match.arg(type)
      
      if(all(class(model) %in% c("rma.mv", "rma")) == FALSE) {stop("Sorry, you need to fit a metafor model of class rma.mv or rma")}
      
      data <- get_data_cont(model, mod_cat)
      
      if(length(mod_cont) == 0){
        # Get confidence intervals
        CI <- get_est(model, mod_cat)
        
        # Get prediction intervals
        PI <- get_pred(model, mod_cat)
        
        model_results <- list(mod_table = cbind(CI, PI[,-1]), data = data)
      }
      
      if(length(mod_cont) >= 1){
        ests <- get_pred_est_cont(model, mod_cat, mod_cont, type)
        model_results <- list(mod_table = ests, data = data)	
      }
      
      class(model_results) <- "orchard"
      
      return(model_results)
      
    }
    
    
# Convert proportions to a logit scale making them normally distributed, which should better satisfy assumptions of SMD and lnCVR
convert_propor <- function(study_ID, obs, male_mean, male_sd, female_mean, female_sd){
    logit_mean_male <- gtools::logit(male_mean) + (male_sd^2 / 2)*((1 / (1-male_mean)^2) - (1 / male_mean^2))
     
      sd_logit_male  <- sqrt(male_sd^2* (((1 / male_mean) + (1 / (1-male_mean)))^2))
      
      logit_mean_female <- gtools::logit(female_mean) + (female_sd^2 / 2)*((1 / (1-female_mean)^2) - (1 / female_mean^2))
      
      sd_logit_female  <- sqrt(female_sd^2* (((1 / female_mean) + (1 / (1-female_mean)))^2))

    return(data.frame(study_ID = study_ID, obs = obs, male_mean_conv = logit_mean_male, male_sd_conv = sd_logit_male, female_mean_conv = logit_mean_female, female_sd_conv = sd_logit_female))
}

## Just a quick check to see things are sensible
# propmale = c(0.20, 0.5, 0.6)
# sdmale = c(0.10, 0.2, 0.4)
# propfemale = c(0.30, 0.6, 0.2)
# sdfemale = c(0.2, 0.2, 0.4)
# convert_propor(prop, sd)

# Convert latency data to a log normal distributed which should normalize and better satisfy assumptions of SMD and lnCVR
convert_latency <- function(study_ID, row_ID, male_mean, male_sd, female_mean, female_sd, method = c("delta", "analyt")){

  if(method == "delta"){
         lnMean <- log(mean) - (sd^2 / (2*mean^2))
      sd_lnMean <- sqrt(sd^2 / mean^2)
  }

  if(method == "analyt"){
         lnMean_male <- log(male_mean) - log(sqrt(1 + (male_sd^2 / male_mean^2)))
      sd_lnMean_male <- sqrt(log(1 + (male_sd^2 / male_mean^2)))
      lnMean_female <- log(female_mean) - log(sqrt(1 + (female_sd^2 / female_mean^2)))
      sd_lnMean_female <- sqrt(log(1 + (female_sd^2 / female_mean^2)))
  }

  return(data.frame(study_ID = study_ID, row_ID = row_ID, male_mean_conv = lnMean_male, male_sd_conv = sd_lnMean_male, female_mean_conv = lnMean_female, female_sd_conv = sd_lnMean_female))
}

# Check that the latency calculation is correct as I have never applied conversions before
# mean_norm <- rnorm(100000, 5, 1)  # Simulate a large N from normal dist

# Calculate mean and sd from this as if it were a sample; should be very precise because of the large N
#mean_normal_scale <- mean(mean_norm)
#  sd_normal_scale <-   sd(mean_norm)

# Now, this is what we should get for mean and sd on a log normal scale
#mean_lognormal_scale <- mean(log(mean_norm))
#  sd_lognormal_scale <-   sd(log(mean_norm))

# Check that our function recovers this, and different methods are nearly the same
#convert_latency(mean_normal_scale, sd_normal_scale, method = "analyt")
#convert_latency(mean_normal_scale, sd_normal_scale, method = "delta")

# These virtually match, so we know the solutions are correct. There will be slight differences given that the simulated data depends on sample size and will vary slight from one sim to the next. The equations are analytical approximations or solutions for what we can expect from repeated simulations.