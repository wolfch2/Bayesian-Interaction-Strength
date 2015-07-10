#########################################################
######################################################### 1. setup (user may have to modify this seciton)
#########################################################

setwd("C:/Users/Chris/Dropbox/stat phd/project_code") # contains folders: "scripts", "data", "models", and "output"
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

#########################################################
######################################################### 3. estimate attack rates!
#########################################################

# load attack rate estimation procedures
methods_files = c("bayes.r", "nonpara_bs.r", "para_bs.r")
sapply(paste0("scripts/", methods_files), source)

# estimate attack rates with bootstrap and Bayesian methods
results = c(list("Nonpara BS" = nonpara_bs(data), "Para BS" = para_bs(data)),
	bayesian_method(data, c=c(0,1/3,1)))

saveRDS(results, "output/results.RDS") # save results for future use (optional)
# results = readRDS("output/results.RDS")

#########################################################
######################################################### 4. figures and graphs
#########################################################

source("scripts/indep_plot.r") # assess relationships between covariates and feeding props.
source("scripts/mle_comp.r") # plots comparing Dirichlet-multinomial post. median to MLE

source("scripts/side_plot.r") # plot comparing point and interval estimates across methods
source("scripts/comp_hists.r") # for Bayes(1/3) method, separate out estimator components

cat("\r") # follow w/ blank line to reset cursor

