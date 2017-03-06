# Each time Bayes model is fit, run it multiple times (on diff. cores), just outputing the first result (so code runs nicely).

###
### Bayesian sampling functions (single prey species)
###

# get handling time model coeffs. using Bayesian regression
get_coeffs = function(y, log.pred.size, log.prey.size, log.temp, len, prey){
	N = length(log.pred.size) # number of observations
	x=cbind(rep(1,N),log.pred.size,log.prey.size,log.temp) # predictors

	data = list(y=y,x=x,N=N,len=len)

	samps = foreach(i=1:250, .packages=c("rjags"), .export=c("nsamples")) %dopar% {
		inits = list(beta=rnorm(4),tau=exp(1),z=runif(N,0,1)) # initial values
		jags.m = jags.model(file = "models/regression.jags", data = data,
			inits = inits, n.chains = 1, n.adapt = 100)
		update(jags.m, n.iter=10000) # burn-in
		params = c("beta") # parameter (vector) to monitor
		samps = coda.samples(jags.m, params, n.iter = nsamples*200, thin=200)
	}
	saveRDS(samps, paste0("output/coeffs_" , prey, ".RDS"))
	
	return(as.matrix(samps[[1]]))
}

# use iid multivariate normal likelihood to get mean and covariance matrix parameters for covariates
get_covar_means = function(log.pred.size, log.prey.size, log.temp, prey){
	y = cbind(log.pred.size, log.prey.size, log.temp)
	N = length(log.pred.size)

	# Priors (relatively non-informative)
	mu0 = as.vector(c(0,0,0))
	S2 = diag(3)/1000
	S3 = diag(3)/10000
	I3 = diag(3)*10
	
	data = list(y=y,N=N,S2=S2,S3=S3,mu0=mu0)

	samps = foreach(i=1:250, .packages=c("rjags","MASS","MCMCpack"), .export=c("nsamples")) %dopar% {
		inits = list(mu = mvrnorm(1,mu0,I3), tau = rwish(3,matrix(diag(3),nrow=3))) # initial values
		jags.m = jags.model(file = "models/multivariate_normal.jags", data = data,
			inits = inits, n.chains = 1, n.adapt = 100)
		update(jags.m, n.iter=10) # burn-in
		params = c("mu", "Sigma") # parameteres of interest
		samps = coda.samples(jags.m, params, n.iter = nsamples, thin=1)
	}
	saveRDS(samps, paste0("output/covars_" , prey, ".RDS"))
	
	samps_mat = as.matrix(samps[[1]])
	return(list(means=samps_mat[,10:12],vars=samps_mat[,1:9]))
}

# Bayes iid gamma model (prob. of zero is handled in calling function with beta prior)
bayes_gamma = function(y, prey){
	N = length(y)
	data = list(y=y,N=N)

	# probably have to export "nsamples" since it's found in the calling environment rather that the function itself..
	# getting error on log.beta w/ very spread out inits (e.g. rnorm(1)): https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/fdbe7c22/?limit=25
	samps = foreach(i=1:250, .packages=c("rjags"), .export=c("nsamples")) %dopar% {
		inits = list(log.alpha=runif(1,0.9,1.1)*log(mean(y)^2/var(y)), log.beta=runif(1,1-1e-2,1+1e-2)*log(mean(y)/var(y)))
		jags.m = jags.model(file = "models/gamma.jags", data = data,
			inits = inits, n.chains = 1, n.adapt = 100)
		update(jags.m, n.iter=1000) # burn-in
		params = c("alpha", "beta")
		samps = coda.samples(jags.m, params, n.iter = nsamples*10, thin=10)
	}
	saveRDS(samps, paste0("output/density_" , prey, ".RDS"))
	
	samps_mat = as.matrix(samps[[1]])
	return(samps_mat[,1]/samps_mat[,2])
}

###
### Bayesian sampling functions (multiple prey species)
###

# get coeffs. for handling time regression model using lab data
get_coeffs_all = function(data){
	# matrix to hold sampled betas:
	coeffs_mat =  array(NA, dim=c(length(data$prey_types),nsamples,4))
	dimnames(coeffs_mat)[[1]] = data$prey_types

	for(prey in data$prey_types){
		if(prey == "Austrolittorina antipodum") # assumed same as Aus. cin.
			next
		htimes_prey = data$htimes[which(data$htimes$Prey == prey),] # relevant handling times
		htimes_prey = na.omit(htimes_prey[,c(4,8,11,26,27,28)]) # remove NAs
		log.pred.size = log(htimes_prey$PredSize)
		log.prey.size = log(htimes_prey$PreySize)
		log.temp = log(htimes_prey$Temperature)
		len = htimes_prey$Max_.days. - htimes_prey$Min_.days. # compute censor interval lengths
		y = htimes_prey$Midpoint_.days. # use window midpoints as responses
		
		if (length(y) == 1){ # use average if only 1 obs.
			coeffs_mat[prey,,] = 0
			coeffs_mat[prey,,1] = y
		} else{
			coeffs_mat[prey,,] = get_coeffs(y, log.pred.size, log.prey.size, log.temp, len, prey)
			if (prey == "Austrolittorina cincta") # Since Aus ant. lacks htime obs. -- see unique(htimes_data$htimes$Prey)
				coeffs_mat["Austrolittorina antipodum",,] = get_coeffs(y, log.pred.size, log.prey.size, log.temp, len, prey)
		}
	}
	return(coeffs_mat)
}

# Sample From Covariate Means and then from Mean handling times
get_htimes = function(data, coeffs_mat){

	mean.h.times = array(NA, dim=c(nsamples, length(data$prey_types)))
	colnames(mean.h.times) = data$prey_types

	for(prey in data$prey_types){
		# get mean/var for covariates
	    	feed_prey = data$feed[which(data$feed$Prey == prey),]
		log.pred.size = log(feed_prey$PredSize)
		log.prey.size = log(feed_prey$PreySize)
		log.temp = log(feed_prey$MeanTemp)		
		out = get_covar_means(log.pred.size, log.prey.size, log.temp, prey) # Sample from mean handling times
		means = out$means
		vars = out$var

		# combine w/ regression coeffs. using post. pred. dist. (and WLLN) -- see paper
		for(i in 1:nsamples){
			betas = coeffs_mat[prey,i,]
			mu = means[i,]
			Sigma = matrix(vars[i,], nrow=3, ncol=3)
			mean.h.times[i,prey] = mean(exp(cbind(1, mvrnorm(1000, mu, Sigma, tol=.001)) %*% betas))
		}
	}
	return(mean.h.times)
}

# multinomial-Dirichlet model (conjugate prior)
get_survey_params_all = function(data, c){
	survey_data = c(data$total_surveyed-sum(data$feeding_status),data$feeding_status)
	samps = rdirichlet(nsamples, survey_data+c)
	return(samps)
}

# sample from mean abundances posterior using ZIG model
get_density_means_all = function(data){
	N_mat = array(NA, dim=c(nsamples, length(data$prey_types)))
	colnames(N_mat) = data$prey_types

	for(prey in data$prey_types){
		y = data$quads[,prey] # abundances for this prey species
		p = rbeta(nsamples, 1 + sum(y==0), 1 + length(y) - sum(y==0)) # sample from post. prob. of a zero
		N_mat[,prey] = (1-p) * bayes_gamma(y[y>0], prey) # treat non-zero parts as iid gamma
	}
	return(N_mat)
}

###
### Bayesian MCMC
###

bayesian_method = function(data, c=1/3){
	coeffs_mat = get_coeffs_all(data) # regression model coeffs.
	mean.h.times = get_htimes(data, coeffs_mat) # combine coeffs. and covariates information

	N_mat = get_density_means_all(data) # abundances

	AR_mats = lapply(c, function(x){
		p_mat = get_survey_params_all(data, x) # feeding proportions (dependent on sparsity prior choice)

		if(x==1/3){ # save components for Bayes(1/3) to plot later
		    	comp_list = list(N_mat=N_mat, mean.h.times=mean.h.times,
				ratio_mat = p_mat[,2:ncol(p_mat)] / matrix(p_mat[,1],nrow=nsamples,ncol=length(data$prey_types),byrow=FALSE))
			saveRDS(comp_list, "output/comp_list.RDS")
		}

		# combine components to obtain samples from attack rates' posterior
		AR_mat = p_mat[,2:ncol(p_mat)] / matrix(p_mat[,1],nrow=nsamples,ncol=length(data$prey_types),byrow=FALSE) / (mean.h.times * N_mat)
		colnames(AR_mat) = data$prey_types
		return(AR_mat)
	})

	names(AR_mats) = paste0("Bayes (c=", round(c,3), ")") # return a list (indexed by c) of attack rate posterior sample matrices
	return(AR_mats)
}

