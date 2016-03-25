# Separate script showing how attack rate estimation procedure can be applied to subsets of the feeding data.

#########################################################
######################################################### 1. setup (user may have to modify this seciton)
#########################################################

setwd("C:/Users/Chris/Dropbox/stat phd/project_code - whelk update") # contains folders: "scripts", "data", "models", and "output"
nsamples = 1000 # number of posterior (and bootstrap) samples to obtain

# use install.packages(...) if any of these packages are not already installed
require(MASS) # need for mvrnorm
require(MCMCpack) # need for rwish
require(mvnmle) # needed to est. multivar norm params.
require(mvtnorm)
require(numDeriv)
require(Matrix)
require(mcmc)
require(rjags) # call jags from R
require(coda) # analyze output in R
require(ggplot2)
require(reshape2)
require(grid)
require(gridExtra)
require(lattice)
require(colorspace)
require(RColorBrewer)
require(plyr)
require(scales)
require(Hmisc)
require(foreach)
require(doSNOW)

#########################################################
######################################################### 2. load and process data
#########################################################

options("stringsAsFactors"=FALSE) # turn off automatic factors on data import

raw_data = list(
	feed = read.csv("data/NZ-FeedingObservations-Select.csv"),
	htimes_reg  = read.csv("data/NZ-HandlingTimes-MRegnCoeff-MeasuredANDMatched.csv",
		skip=3),
	htimes = read.csv("data/NZ-HandlingTimes.csv"),
	temps = read.csv("data/NZ-Temps-Monthly.csv"),
	quads = read.csv("data/NZ-SpeciesDensities-Select.csv")
)

# prepare feeding survey, handling time, and abundance data for analysis
source("scripts/process_data.r")
data = process_data(raw_data)

quantile_table = table(data$feed$Prey, cut(data$feed$PredSize, quantile(data$feed$PredSize, seq(0,1,by=1/8))))
quantile_table = rbind(quantile_table, `Total Surveyed`=apply(quantile_table,2,sum))
tmp = latex(quantile_table, file = "output/quantile_table.tex", table.env=FALSE, booktabs=TRUE)

# split feeding survey data by predator size into eight groups (of roughly equal size)
feed_split = split(data$feed, cut(data$feed$PredSize, quantile(data$feed$PredSize, seq(0,1,by=1/8))))

# load attack rate estimation procedures (just using Bayes here)
methods_files = c("bayes.r", "nonpara_bs.r", "para_bs.r")
sapply(paste0("scripts/", methods_files), source)

registerDoSNOW(makeCluster(8, type = "SOCK"))
getDoParWorkers()

pack_list = c("MASS","MCMCpack","mvnmle","mvtnorm","numDeriv","Matrix","mcmc","coda","rjags")
results_list = foreach(i = 1:length(feed_split), .combine=c, .packages=pack_list) %dopar% {
       data$feed = feed_split[[i]]
       results = bayesian_method(data) # apply Bayesian method to each subset separately
}

names(results_list) = sapply(feed_split,function(x) median(x$PredSize)) # median pred sizes in each group

saveRDS(results_list, "output/results_list_8classes.RDS") # save results for future use (optional)
# results_list = readRDS("output/results_list_8classes.RDS")

#########################################################
######################################################### plot
#########################################################

data = do.call("rbind",sapply(names(results_list),function(x) data.frame(Size=x,results_list[[x]][,c("Chamaesipho columna","Xenostrobus pulex")]), simplify=FALSE))
data = melt(data,id.vars="Size")
data$value = log(data$value,base=10)
data$variable = gsub("\\.", " ", data$variable)
data$Mass = 1.214*10^-4 * as.numeric(as.character(data$Size))^3.210

dev.new(width=7,height=6)
ggplot(data, aes(x=Mass,y=value,group=Mass)) +
	geom_violin(fill="gray") +
	facet_wrap(~ variable, scales="free_y") +
	theme_bw()  +
	stat_summary(aes(group=Mass), fun.y=median, geom="point") +
	scale_alpha(range=c(0,1), guide=FALSE) +
	xlab("Median predator mass (g)") +
	theme(strip.text = element_text(size = 10, face = "italic", angle = 0)) +
	scale_y_continuous(labels=trans_format("I", math_format(10^.x))) +
	scale_x_log10(breaks=c(.1,.5,1)) +
	expand_limits(x=.1) +
	annotation_logticks(sides="b") +
	ylab(bquote("Attack Rate (" ~ xi[i] ~ ")"))

savePlot(filename = "output/size_class.pdf", type = "pdf")
savePlot(filename = "output/size_class.png", type = "png")
dev.off()

