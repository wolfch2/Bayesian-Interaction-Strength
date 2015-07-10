# plot average covariates versus feeding proportions to look for dependence

feed = data$feed
prey_types = data$prey_types

covar_means_mat =  array(NA, dim=c(length(prey_types),nsamples,4))
dimnames(covar_means_mat)[[1]] = prey_types

covariates = c("PredSize", "PreySize", "MeanTemp")
covariate_names_short = c("predsize", "preysize", "temp")

# http://r.789695.n4.nabble.com/Degree-symbol-td3415813.html
covariate_names = c(expression(paste("log"[10], " Predator Size (mm)")), 
	expression(paste("log"[10], " Prey Size (mm)")), 
	expression(paste("log"[10], " Temperature (",degree,"C)")))

df = NULL

# combine avg. (log) covar and feeding prop. info. for species with at least 3 obs.
for(covar in 1:3){
	for(prey in sort(prey_types)){
		prop = NULL
		covar_vec = NULL
		for(i in unique(feed$CensusID)){
			obs = intersect(which(feed$Prey == prey), which(feed$CensusID == i))
			if(length(obs) == 0)
				next
			prop = c(prop, length(obs)/sum(feed$CensusID == i))
			covar_vec = c(covar_vec, mean(log(unlist(     feed[obs,][covariates[covar]]     ),base=10)))
		}
		if(length(prop) > 3){
			df = rbind(df, cbind.data.frame(prey,covar,covar_vec,prop))
		}
	}
}
df$covar = factor(df$covar, levels=1:3, labels=covariate_names)
df$prey = gsub(" ", "\n", df$prey)

# handle label formatting
plot_labeller <- function(variable,value){
	if (variable=="prey"){
		return(value)
	} else{
		return(covariate_names[value])
	}
}

ggplot(as.data.frame(df), aes(x=covar_vec,y=prop)) +
	geom_point() +
	facet_grid(prey~covar, scales="free_x", labeller=plot_labeller) +
	ylab("Predators Feeding") +
	xlab("") +
	theme_bw() +
	theme(panel.margin = unit(1, "lines"),
		strip.text.y = element_text(size = 10, face = "italic", angle = 0)) +
	scale_y_continuous(labels = percent)

ggsave("output/indep.pdf")
dev.off()
