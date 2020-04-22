
#model=model4.RR
#mod_cat = "t_magnitude"
#mod_cont="mean_t"
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

#test <- mod_results(model5_RR, mod_cat = "order", mod_cont = c("mean_t", "delta_t", "Q10"))
#print(test)

#orchard_plot(test, mod = "order", xlab = "Order", angle=45) 