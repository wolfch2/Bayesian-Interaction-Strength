######################################################## Parametric bootstrap

# Strategy -- estimate all the parameters (mu's, beta's, ...) needed to generate
# new data sets.  Then, need only gen. new data, compute AR ests., take quantiles.

##################### Get Handling Times Data, Est Coeffs.

get_coeffs_mat = function(data){
	# matrix to hold estimated betas:
	coeffs_mat =  array(NA, dim=c(length(data$prey_types),5))
	dimnames(coeffs_mat)[[1]] = data$prey_types

	for(prey in data$prey_types){
	    	# comments on data processing are in Bayes.r
		if(prey == "Austrolittorina antipodum") # assumed same as Aus. cin.
			next
		htimes_prey = data$htimes[which(data$htimes$Prey == prey),]
		htimes_prey = na.omit(htimes_prey[,c(4,8,11,26,27,28)])
		log.pred.size = log(htimes_prey$PredSize)
		log.prey.size = log(htimes_prey$PreySize)
		log.temp = log(htimes_prey$Temperature)
		y = log(htimes_prey$Midpoint_.days.) # interval censoring not accounted for here, response log transformed (similar to NPBS)
		
		if (length(y) == 1){ # use average if only 1 obs. (i.e. sd = 0 -- wrong but no information)
			coeffs_mat[prey,] = 0
			coeffs_mat[prey,1] = y
		} else{
		    	# estimate parameters with maximum likelihood (lm)
			coeffs_mat[prey,1:4] = lm(y~log.pred.size+log.prey.size+log.temp)$coefficients
			coeffs_mat[prey,5] = summary(lm(y~log.pred.size+log.prey.size+log.temp))$sigma
			if (prey == "Austrolittorina cincta"){ # use this model for Aus. ant. as well
				coeffs_mat["Austrolittorina antipodum",1:4] = lm(y~log.pred.size+log.prey.size+log.temp)$coefficients
				coeffs_mat["Austrolittorina antipodum",5] = summary(lm(y~log.pred.size+log.prey.size+log.temp))$sigma
			}
		}
	}
	return(coeffs_mat)
}

##################### get abundance ZIG parameters

# negative log likelihood for iid gamma
neg.log.lik = function(theta, y){
	return(-sum(log(dgamma(y, theta[1], theta[2]))))
}

# iid gamma maximum likelihood estimate
gamma_mle = function(y){
	return(nlminb(start=c(sqrt(mean(y)), sqrt(mean(y))), objective=neg.log.lik, y=y)$par)
}

# Resample prey species densities (abundances) using parametric bootstrap.
#
# cols. of N_mta are p (MLE for prob. of zero under binom. likelihood -- just sample prop.)
# followed by gamma dist. MLE ests. of alpha and beta.
# rows are prey species.
resample_densities = function(data, N_mat){
	density = matrix(NA, nrow=nrow(data$quads), ncol=length(data$prey_types))
	colnames(density) = data$prey_types
	for(prey in data$prey_types)
		density[,prey] = rbinom(nrow(data$quads), 1, 1-N_mat[prey,1]) * rgamma(nrow(data$quads), N_mat[prey,2], N_mat[prey,3])
	return(density)
}

# estimate multivariate normal paramters for (log) covariates dist.
get_covariates_dist_ests = function(data){
	covar_means_mat = array(NA, dim=c(length(data$prey_types),3))
	dimnames(covar_means_mat)[[1]] = data$prey_types

	covar_covar_mat = array(0, dim=c(length(data$prey_types),3,3)) # covariates covariance array
	dimnames(covar_covar_mat)[[1]] = data$prey_types

	for(prey in data$prey_types){
		feed_prey = data$feed[which(data$feed$Prey == prey),]
		log.pred.size = log(feed_prey$PredSize)
		log.prey.size = log(feed_prey$PreySize)
		log.temp = log(feed_prey$MeanTemp)

		# MLE est. for mean is vector of sample proportions
		covar_means_mat[prey,] = c(mean(log.pred.size), mean(log.prey.size), mean(log.temp))

		# keep sig hat = 0 if length == 1 (can't do much else)
		# if fewer than 10 obs., just estimate diagonal (assume components independent)
		if(1 < length(log.temp) && length(log.temp) < 10){ # hard to est. sig hat here
			covar_covar_mat[prey,1,1] = var(log.pred.size)
			covar_covar_mat[prey,2,2] = var(log.prey.size)
			covar_covar_mat[prey,3,3] = var(log.temp)
		} else if(length(log.temp) >= 10){
			covar_covar_mat[prey,,] = mlest(cbind(log.pred.size, log.prey.size, log.temp))$sigmahat
		}
	}
	return(list(covar_means_mat=covar_means_mat, covar_covar_mat=covar_covar_mat))
}

para_bs = function(data){
    	# matrix to hold attack rate estimates
	AR_mat =  array(NA, dim=c(length(data$prey_types), nsamples))
	dimnames(AR_mat)[[1]] = data$prey_types

	# estimate feeding survey parameters (using sample proportions)
	p.vec = data$feeding_status/data$total_surveyed

	# estimate abundance parameters
	N_mat = array(NA, dim=c(length(data$prey_types), 3)) # parameter order: p, alpha, beta (rate)
	dimnames(N_mat)[[1]] = data$prey_types
	for(prey in data$prey_types){
		y = as.numeric(as.matrix(data$quads[,prey]))
		N_mat[prey, 1] = mean(y < 10^-5) # use sample prop. of zeros as est. for prob. of zero
		if(sum(y < 10^-5) > 0)
			y = y[-which(y < 10^-5)]
		ml = gamma_mle(y/10^6) # gamma MLE for remaining part
		ml[2] = ml[2]/10^6
		N_mat[prey, 2:3] = ml 
	}

	# estimate reg. coeff. and covariate dist. parameters 
	covar_dist_ests = get_covariates_dist_ests(data)
	coeffs_mat = get_coeffs_mat(data)

	# resample new datasets and compute attack rate estimates
	for(i in 1:nsamples){
		density = apply(resample_densities(data, N_mat), 2, mean) # resample abundances
		bs.survey = rmultinom(1, data$total_surveyed, c(1-sum(p.vec), p.vec))[-1] # resample feeding survey using sample props.
		A = bs.survey / data$total_surveyed # sample proprotions for resampled survey (to compute attack rate estimates)
		names(A) = data$prey_types

		# now simulate avg. handling time for each prey species
		for(prey in data$prey_types){
			if(A[prey] < 10^-10){
				AR_mat[prey,i] = 0 # AR est is 0 if none of this species is observed in resampled survey
				next
			}
			# resample covariates
			covars = mvrnorm(n=data$total_surveyed, mu=covar_dist_ests$covar_means_mat[prey,], Sigma=covar_dist_ests$covar_covar_mat[prey,,])
			coeffs = matrix(coeffs_mat[prey,1:4], nrow=data$total_surveyed, ncol=4, byrow=TRUE)
			coeffs * cbind(1, covars)
			meanhtime = mean(exp(apply(coeffs * cbind(1, covars),1,sum))) # mean handling time based on resampled covariates
			AR_mat[prey,i] = A[prey] / ( (1-sum(A)) * meanhtime * density[prey] ) # resampled attack rates
		}
	}
	return(t(AR_mat))
}

