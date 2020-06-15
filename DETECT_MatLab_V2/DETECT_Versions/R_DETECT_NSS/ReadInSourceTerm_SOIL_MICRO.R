print("Reading in params")
params <- read.csv(file="Inputs2/ST_params.csv", header=FALSE)
print("###################")

print("Reading in params_lag")
params_lag <- read.csv(file="Inputs2/ST_params_lag.csv", header=FALSE)
print("###################")

print("Reading in Xrel")
Xrel <- read.csv(file="Inputs2/ST_Xrel.csv", header=FALSE)
print("###################")

print("Reading in SWC")
SWC <- read.csv(file="Inputs2/ST_SWC.csv", header=FALSE)
print("###################")

print("Reading in SoilT")
SoilT <- read.csv(file="Inputs2/ST_SoilT.csv", header=FALSE)
print("###################")

print("Reading in SWCant_Rm")
SWCant_Rm <- read.csv(file="Inputs2/ST_SWCant_Rm.csv", header=FALSE)
print("###################")

print("Reading in SWCant_Rr")
SWCant_Rr <- read.csv(file="Inputs2/ST_SWCant_Rr.csv", header=FALSE)
print("###################")

print("Reading in SoilTant")
SoilTant <- read.csv(file="Inputs2/ST_SoilTant.csv", header=FALSE)
print("###################")

print("Reading in Gness")
Gness <- read.csv(file="Inputs2/ST_Gness.csv", header=FALSE)
print("###################")

###            R = [zeros(Nt,1) Rall.R];
print("Reading in R")
R <- read.csv(file="Inputs3/R_R.csv", header=FALSE)
print("###################")

###            Rr= [zeros(Nt,1) Rall.Rr];
print("Reading in Rr")
Rr <- read.csv(file="Inputs3/R_Rr.csv", header=FALSE)
print("###################")

###            Rm= [zeros(Nt,1) Rall.Rm];   
print("Reading in Rm")
Rm <- read.csv(file="Inputs3/R_Rm.csv", header=FALSE)
print("###################")

