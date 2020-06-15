#read_file("Dataset2b_SI.txt")
remove(list = ls())
setwd(here:::here())
library(R2OpenBUGS)
library(coda)
library(readr)
library(R.matlab) # MAY BE NOT USED IN FUTURE/ CHECK IF REMOVAL NEEDED

# Read all processed data in
source("DETECT_R_BETA_V1_0_PLUS_BAYES.R")
# TEMP alternative read in
# Nt <- 732
# Nz <- 100
# Stype <- as.matrix(readMat("Stype_mat.mat"))
#Stype <- read.csv(file="Inputs4/Stype_data.csv", header=FALSE)
# Dgs_t <- as.matrix(read.csv(file="Inputs4/Dgs_t_data.csv", header=FALSE))
# class(Dgs_t)
# mode(Dgs_t)
# catm <- as.matrix(read.csv(file="Inputs4/catm_data.csv", header=FALSE))
# class(catm)
# mode(catm)
# 
# attributes(Stype)
# List of data feeded to the model
datlist <- list(Nt=Nt,
                Nz = Nz,
                Stype = Stype, 
                Dgs_t = Dgs_t,
                catm = upbc)


# List of inital values
#modelInits <- list(list(),list(),list())

modelInits <- list(list(tau_cctype = 0.1, cctype = 0.1),list(tau_cctype = 0, cctype = 0),list(tau_cctype = 0.01, cctype = 0.01))

# Parameters
modelParameters<- c("tau_cctype", "cctype")


# Number of model iterations
modelIterations <- 5000

# Number of chains
modelChains <- 3

# Number of iterations which considered burnin/warmup
modelBurnin <- 1000

model_run  <-  bugs(data=datlist,
                    inits=modelInits,
                    parameters.to.save=modelParameters,
                    n.iter=modelIterations,
                    n.chains=modelChains, n.burnin=modelBurnin, n.thin=1,
                    model.file="RAW_DETECT_OPENBUGS.R",
                    codaPkg = T,debug=T)
