#########################################################
######################################################### 1. setup (user may have to modify this seciton)
#########################################################

setwd("/home/stats/wolfch/project_code") # contains folders: "scripts", "data", "models", and "output"
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
require(scales)
require(ggplot2)
require(reshape2)
require(grid)
require(gridExtra)
require(lattice)
require(colorspace)
require(RColorBrewer)
require(plyr)
require(foreach)
require(doSNOW)

cl = makeCluster(25, type = "SOCK", outfile="/home/stats/wolfch/project_code/worker_log.txt")
registerDoSNOW(cl)

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

# load attack rate estimation procedure
source("scripts/bayes.r")

# estimate attack rates with bootstrap and Bayesian methods
results = bayesian_method(data, c=1/3)

# saveRDS(results, "output/results.RDS") # save results for future use (optional)
# results = readRDS("output/results.RDS")

stopCluster(cl)