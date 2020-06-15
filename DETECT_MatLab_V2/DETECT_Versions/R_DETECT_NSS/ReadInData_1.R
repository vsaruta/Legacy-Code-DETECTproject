print("Reading in data required for DETECT model")

print("Reading in Atmoshperic pressure")
AtmPress <- read.csv(file="Inputs/DD_AtmPress.csv", header=FALSE)
#View(AtmPress)
print("###################")

print("Reading in BC")
BC <- read.csv(file="Inputs/DD_BC.csv", header=FALSE)
#View(BC)
print("###################")

print("Reading in d13Cm")
d13Cm <- read.csv(file="Inputs/DD_d13Cm.csv", header=FALSE)
#View(d13Cm)
print("###################")

print("Reading in d13Cr")
d13Cr <- read.csv(file="Inputs/DD_d13Cr.csv", header=FALSE)
#View(DD_d13Cr)
print("###################")

print("Reading in dx")
dx <- read.csv(file="Inputs/DD_dx.csv", header=FALSE)
#View(dx)
print("###################")

print("Reading in Gness")
Gness <- read.csv(file="Inputs/DD_Gness.csv", header=FALSE)
#View(Gness)
print("###################")

print("Reading in ic")
ic <- read.csv(file="Inputs/DD_ic.csv", header=FALSE)
#View(ic)
print("###################")

print("Reading in Nt_run")
Nt_run <- read.csv(file="Inputs/DD_Nt_run.csv", header=FALSE)
#View(Nt_run)
print("###################")

print("Reading in params")
params <- read.csv(file="Inputs/DD_params.csv", header=FALSE)
#View(params)
print("###################")

print("Reading in params_lag")
params_lag <- read.csv(file="Inputs/DD_params_lag.csv", header=FALSE)
#View(params_lag)
print("###################")

print("Reading in params_swp")
params_swp <- read.csv(file="Inputs/DD_params_swp.csv", header=FALSE)
#View(params_swp)
print("###################")

print("Reading in SoilT")
SoilT <- read.csv(file="Inputs/DD_SoilT.csv", header=FALSE)
#View(SoilT)
print("###################")

print("Reading in SoilTant")
SoilTant <- read.csv(file="Inputs/DD_SoilTant.csv", header=FALSE)
#View(SoilTant)
print("###################")

print("Reading in SWC")
SWC <- read.csv(file="Inputs/DD_SWC.csv", header=FALSE)
#View(SWC)
print("###################")

print("Reading in SWCant_Rm")
SWCant_Rm <- read.csv(file="Inputs/DD_SWCant_Rm.csv", header=FALSE)
#View(SWCant_Rm)
print("###################")

print("Reading in SWCant_Rr")
SWCant_Rr <- read.csv(file="Inputs/DD_SWCant_Rr.csv", header=FALSE)
#View(SWCant_Rr)
print("###################")

print("Reading in Xrel")
Xrel <- read.csv(file="Inputs/DD_Xrel.csv", header=FALSE)
#View(Xrel)
print("###################")

print("###################")
print("###################")