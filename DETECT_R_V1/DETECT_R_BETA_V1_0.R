# % Based on soil CO2 production model developed by Fang & Moncrieff (1999)
# % Agricultural and Forest Meteorology 95:225-236. But, simplify model to
# % only deal with gas phase CO2 as is done in Hui & Luo (2004) Global
# % Biogeochemical Cycles Vol 18:GB4029.
# % Translated in to "R" by Volodymyr Saruta (vsaruta@gmail.com)
# %
# % INPUTS:
# % Nt_sim = no. of time-steps to run the model for.  
# % params, Xrel, SWC, SoilT, SWCant, SoilTant, Gness are used in SourceTerm.m.
# % Xrel stores the distr. by depth of soil CO2, microb. C & root C (mgC/cm-3).
# % SWC(t,z) = soil water content (m3/m3) at time t, depth z
# % SoilT(t,z) = soil temp (C) at time t, depth z
# % SWCant(t,z,lag) = soil water content (m3/m3) at time t - lag, depth z
# % SoilTant(t,z) = soil temp (C) at time t - lag, depth z
# % Gness(t,1) = Vegetation greeness at time t. 
# % AtmPress(t,1) = atmostpheric pressure at the surface of PHACE.
# % phig is the air filled porosity (m3/m3); phig100 is phig at -100cm H20.
# % dx=[dt, dz]
# % BC = [cppm, cdeep] (boundary conditions), where cppm = atmospheric [CO2]
# % (or [CO2] at the soil surface) (ppm); cdeep = soil [CO2] at the bottom
# % soil node (mg CO2 m-3).
# % cinit(1,z) = initial conditions (soil [CO2], mg CO2 m-3) at time "zero"
# % d13Cr = d13C of root-respired CO2, dim = 1xNz
# % d13Cm = d13C of mirobial-respired CO2, dim = 1xNz
# %
# % OUTPUTS:
# % out.CO2type =  is a matrix of [CO2] for every time (row), at every depth
# % (col.) and CO2 type (12C-micr, 13C-micr, 12C-roots, 13C-roots).
# % out.CO2 = the total [CO2], i.e. we sum across all four types.
# % out.CO2fluxtype is a matrix of soil surface CO2 flux for every type (row)
# % and time (col.).
# % out.CO2flux = the total soil surface CO2 flux, i.e. we sum across all four types.
# % out.d13C = total d13C of CO2.
# % out.roots = soil surface CO2 flux from roots.

print("###################")
print("Initializing script")
print("###################")


print("Clearing enviroment")
rm(list=ls());print("Enviroment is cleared");
remove(list = ls())
print("###################")

print("Setting up working directory")
setwd(here:::here());print("Working directory is changed to:"); getwd()
print("###################")



print("Fetching model run properties")
running_antecedent_vers <- 0;  #% Including (=1) or not (=0) the antecedent variables in the DETECT model.
runmodel <- 1;                 #% Whether to run the model (=1) or not(=0).
Nt_run <- 183*4;               #% The number of time steps the user wishes to run the DETECT model for. For 
#% Ryan et al. (in review), max(Nt_run)=183*4, where 183 is the no. of days 
#% model is run for and 4 is the no. of 6 hourly time-steps each day.
Nt_beg<- 92;                  #% First DoY to use for computer experiment.
Nt_end<- 274;                 #% Final DoY to use for computer experiment.
writefile<- 1;                #% Whether to write the output to mat and csv files (=1) or not (=0).
print("###################")

print("Adding required libraries/packages")
library(rmatio)
library(R.matlab)
library(pracma)
library(tidyverse)
library(imager)
print("###################")

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
# SoilTant <- read.csv(file="Inputs/DD_SoilTant.csv", header=FALSE)
SoilTant <-  sample(0,292800,replace = TRUE)
dim(SoilTant) <- c(732,100,4)
SoilTant_Array <- read.csv(file="Inputs/DD_SoilTant.csv", header=FALSE)
SoilTant_1 <- SoilTant_Array[,1:100]
SoilTant_2 <- SoilTant_Array[,101:200]
SoilTant_3 <- SoilTant_Array[,201:300]
SoilTant_4 <- SoilTant_Array[,301:400]
# SoilTant_1 <- read.csv(file="Inputs/DD_SoilTant.csv", header=FALSE)[,1:100]
# SoilTant_2 <- read.csv(file="Inputs/DD_SoilTant.csv", header=FALSE)[,101:200]
# SoilTant_3 <- read.csv(file="Inputs/DD_SoilTant.csv", header=FALSE)[,201:300]
# SoilTant_4 <- read.csv(file="Inputs/DD_SoilTant.csv", header=FALSE)[,301:400]
for (i in 1:732)
{
  for (j in 1:100)
  {
    SoilTant[i,j,1] <- SoilTant_1[i,j]
    SoilTant[i,j,2] <- SoilTant_2[i,j]
    SoilTant[i,j,3] <- SoilTant_3[i,j]
    SoilTant[i,j,4] <- SoilTant_4[i,j]
  }
}
print("###################")


print("Reading in SWC")
SWC <- read.csv(file="Inputs/DD_SWC.csv", header=FALSE)
print("###################")

print("Reading in SWCant_Rm")
SWCant_Rm <-  sample(0,292800,replace = TRUE)
dim(SWCant_Rm) <- c(732,100,4)
SWCant_Rm_Array <- read.csv(file="Inputs/DD_SWCant_Rm.csv", header=FALSE)
SWCant_Rm_1 <- SWCant_Rm_Array[,1:100]
SWCant_Rm_2 <- SWCant_Rm_Array[,101:200]
SWCant_Rm_3 <- SWCant_Rm_Array[,201:300]
SWCant_Rm_4 <- SWCant_Rm_Array[,301:400]
# SWCant_Rm_1 <- read.csv(file="Inputs/DD_SWCant_Rm.csv", header=FALSE)[,1:100]
# SWCant_Rm_2 <- read.csv(file="Inputs/DD_SWCant_Rm.csv", header=FALSE)[,101:200]
# SWCant_Rm_3 <- read.csv(file="Inputs/DD_SWCant_Rm.csv", header=FALSE)[,201:300]
# SWCant_Rm_4 <- read.csv(file="Inputs/DD_SWCant_Rm.csv", header=FALSE)[,301:400]
for (i in 1:732)
{
  for (j in 1:100)
  {
    SWCant_Rm[i,j,1] <- SWCant_Rm_1[i,j]
    SWCant_Rm[i,j,2] <- SWCant_Rm_2[i,j]
    SWCant_Rm[i,j,3] <- SWCant_Rm_3[i,j]
    SWCant_Rm[i,j,4] <- SWCant_Rm_4[i,j]
  }
}
print("###################")

print("Reading in SWCant_Rr")
SWCant_Rr <-  sample(0,292800,replace = TRUE)
dim(SWCant_Rr) <- c(732,100,4)

SWCant_Rr_Array <- read.csv(file="Inputs/DD_SWCant_Rr.csv", header=FALSE)
SWCant_Rr_1 <- SWCant_Rr_Array[,1:100]
SWCant_Rr_2 <- SWCant_Rr_Array[,101:200]
SWCant_Rr_3 <- SWCant_Rr_Array[,201:300]
SWCant_Rr_4 <- SWCant_Rr_Array[,301:400]
# SWCant_Rr_1 <- read.csv(file="Inputs/DD_SWCant_Rr.csv", header=FALSE)[,1:100]
# SWCant_Rr_2 <- read.csv(file="Inputs/DD_SWCant_Rr.csv", header=FALSE)[,101:200]
# SWCant_Rr_3 <- read.csv(file="Inputs/DD_SWCant_Rr.csv", header=FALSE)[,201:300]
# SWCant_Rr_4 <- read.csv(file="Inputs/DD_SWCant_Rr.csv", header=FALSE)[,301:400]
for (i in 1:732)
{
  for (j in 1:100)
  {
    SWCant_Rr[i,j,1] <- SWCant_Rr_1[i,j]
    SWCant_Rr[i,j,2] <- SWCant_Rr_2[i,j]
    SWCant_Rr[i,j,3] <- SWCant_Rr_3[i,j]
    SWCant_Rr[i,j,4] <- SWCant_Rr_4[i,j]
  }
}
SWCant_Rr[732,1,1:4]
print("###################")

print("Reading in Xrel")
Xrel <- read.csv(file="Inputs/DD_Xrel.csv", header=FALSE)
print("###################")


print("Running preparation code")



print("Defining constants")
# % Rstar = sum(Rrel) = total root biomass (mgC) mass in a one cm2 column of soil 1m deep. 
Rstar <- params_1 <- 111.54;#  % Rstar = total root biomass (mgC) mass in a one cm2 column of soil 1m deep.

# % RrBase, mass-specific root respiration base rate at a tau degC (mgC cm-3 hr-1).
RrBase <- params_2 <- 0.00006;#  %  RrBase, mass-specific root respiration base rate at a tau degC (mgC cm-3 hr-1).

# % The main effect of SWC on the log-scaled Rrbase; 
alpha1_Rr <- params_3 <- 11.65;  #  % alpha1 for Rr; SWC param.

# % The main effect of SWCant on the log-scaled Rrbase; 
alpha2_Rr <- params_4 <- 20.72; #  % alpha2 for Rr; The main effect of SWCant on the log-scaled Rrbase

# % The main effect of SWC*SWCant on the log-scaled Rrbase; 
alpha3_Rr <- params_5 <- -164.16;#  %  alpha3 for Rr; The main effect of SWC*SWCant on the log-scaled Rrbase

#%Parameters used by the Rm submodel.

# % Cstar in M-M model = total soil C (mgC) mass in a one cm2 column of soil 1m deep.
Cstar <- params_6 <- 711.6094;#  %  Cstar = total soil C (mgC) mass in a one cm2 column of soil 1m deep. 

# % Mstar in M-M model = total MBC (mgC) mass in a one cm2 column of soil 1m deep.
Mstar <- params_7 <- 12.2568;#  % Mstar = total MBC (mgC) mass in a one cm2 column of soil 1m deep. 

# % base respiration rate for microbial respiration function at tau deg. C (mgC cm-3 hr-1).
RmBase <- params_8 <- 0.0015;  #  %  RmBase; (mg C cm-3 hr-1); based on value used in Davidson et al. (2012)

# % The main effect of SWC on the log-scaled Rmbase; 
alpha1_Rm <- params_9 <- 14.05;   #  %  alpha1 (mg C cm-3 hr-1); based on value used in Ryan et al. (2015) 

# % The main effect of SWCant on the log-scaled Rmbase; 
alpha2_Rm <- params_10 <- 11.05;   #  %  alpha2; The main effect of SWCant on the log-scaled Rmbase;  

# % The main effect of SWC*SWCant on the log-scaled Rmbase; 
alpha3_Rm <- params_11 <- -87.55;  #  %  alpha3; The main effect of SWC*SWCant on the log-scaled Rmbase; 

# % Michaelis constants (?) from DAMM paper
Km <- params_12 <- 0.00001;#  %  Km

# % Carbon Use Efficiency (mgC/mgC)
CUE <- params_13 <- 0.8;    #  %  CUE; Carbon Use Efficiency (mgC/mgC).  Taken from Tucker et al. (2012) 

# % Fract. of substrate that Sx that is soluble (called 'p' in DAMM paper). Dimesionless.
solfrac <- params_14 <- 0.0004; #  %  solfrac; Fraction of substrate that Sx that is soluble (called 'p' in DAMM paper). 

# % Diffusivity of substrate Sx that is liquid.  Dimensionless.
Dliq <- params_15 <- 3.17;   #  %  Dliq; Diffusivity of substrate Sx that is liquid.  Dimensionless.

# % Eo term used in Rm submodel (Kelvin)
Eo_Rm <- params_16 <- 324.6;   #  %  Eo_Rm; Eo term used in Rm submodel (Kelvin). Use value from Ryan et al. (2015)

# % Eo term used in Rr submodel (Kelvin)
Eo_Rr <- params_17 <- 324.6;  #  %  Eo_Rr; Eo term used in Rr submodel (Kelvin); Use value from Ryan et al. (2015)

# % Temp. parameter used in the Lloyd & Taylor temp. response type function for Rm (Kelvin)
To <- params_18 <- 227.5;  #  %  To; Temp. parameter used in the Lloyd & Taylor temp. (Kelvin)

# % The main effect of SoilTant on the log-scaled Rmbase Rrbase;
alpha4 <- params_19 <- -4.67;   #  %  alpha4; The main effect of SoilTant on the log-scaled Rmbase;

params_20  <- 10;     #  %  tau; Base temperature in Kelvin.  We use 10 degC.
tau <- params_20+273.15;  # % Base temperature in Kelvin.

if (running_antecedent_vers==0)
{
  alpha2_Rr <- params_4 <- 0; #  % alpha2 for Rr; The main effect of SWCant on the log-scaled Rrbase
  
  # % The main effect of SWC*SWCant on the log-scaled Rrbase; 
  alpha3_Rr <- params_5 <- 0;#  %  alpha3 for Rr; The main effect of SWC*SWCant on the log-scaled Rrbase
  
  #%Parameters used by the Rm submodel.
  # % The main effect of SWCant on the log-scaled Rmbase; 
  alpha2_Rm <- params_10 <- 0;   #  %  alpha2; The main effect of SWCant on the log-scaled Rmbase;  
  
  # % The main effect of SWC*SWCant on the log-scaled Rmbase; 
  alpha3_Rm <- params_11 <- 0;  #  %  alpha3; The main effect of SWC*SWCant on the log-scaled Rmbase; 
  
  # % The main effect of SoilTant on the log-scaled Rmbase Rrbase;
  alpha4 <- params_19 <- 0;   #  %  alpha4; The main effect of SoilTant on the log-scaled Rmbase;
}


# %Lag parameters.  The first row describes the importance of each of hte past
# %4 days on Rm in Source term submodel.  The 2nd row describes the importance 
# %of each of hte past 4 weeks on Rr in Source term submodel.
params_lag <- matrix(data = c(0.75, 0.25, 0, 0,
                              0.2, 0.6, 0.2, 0,
                              0.25, 0.25, 0.25, 0.25), nrow = 3, ncol = 4, byrow = TRUE)

# #%Lag parameters
lagSWC_Rm <- params_lag[1, ];
lagSWC_Rr <- params_lag[2, ];
lagSoilT <- params_lag[3, ];    

# #%The vectors of the relative depth-dependent distributions of microbial biomass
# #%carbon (Mrel), total root biomass (Rrel) and soil organic carbon (Crel).
Crel <- Xrel[1, ];
Mrel <- Xrel[2, ];
Rrel <- Xrel[3, ];

# % Time step (dt) and depth increment (dz)
dt <- dx[1,1]; # % hours (hr)
dz <- dx[1,2]; # % meters (m)

# % Number of time points (Nt), including initial time; and number of soil
# % nodes (Nz), including surface layer (bottom boundary node is Nz+1)
Nt <- nrow(SWC); Nz = ncol(SWC);

# %Create the phig dataset to be used when we compute the diffusion matrix 
# %phig (air filled porosity) = phiT-SWC, where phiT (total soil porosity) = 1-(BD/PD)
# %The value of SWC100 are estimated from fitting the water release curve function to site data of SWP.
BD <- params_swp[1,];
PD <- 2.52; 
SWC100 <- params_swp[3,];
phiT=as.numeric(1-(BD/PD));  # %air-filled porosity
phig=pmin(1,pmax(0,repmat(phiT,Nt,1) - SWC));
phig100=repmat(as.numeric(phiT-SWC100),Nt,1);

# %Predict Diffusion at every depth & time using Moldrup's function. # %Dgs0st = Diffusion 
# %coefficient at standard atmosphere of temperature. T0=273.15K and pressure P0=101.3kPa.
###    Dg0st=repmat(0.0000139,Nt,Nz);  
Dg0st <- repmat(0.0000139,Nt,Nz);

T0 <- 273.15;
b <- params_swp[2,];
P0 <- repmat(101.3,Nt,Nz); # %standard atmospheric pressure
T_T <- SoilT + 273.15;
Atm_Press_as_numeric_unlist <- as.numeric(unlist(AtmPress, use.names=FALSE))
Atm_Press <- repmat(Atm_Press_as_numeric_unlist,1,1);
P <- repmat(t(Atm_Press),1,Nz);
logDg0 <- log(Dg0st)+(1.75*log(T_T/T0))+log(P0/P);
b <- t(as.numeric(unlist(b, use.names=FALSE)))
logDgs <- logDg0+log((2*phig100^3)+(0.04*phig100))+((2+(3/repmat(b,Nt,1)))*log(phig/phig100));

Dgs <- exp(logDgs)*3600;
Dgs_13C <- Dgs/1.0044;
Dgs_12C <- Dgs;

# %Upper boundary condition.  Note that we don't require a lower boundary conditions 
# %because we assume dC/dz = 0 at the lower boundary (i.e. at 1m depth)  
catm <- BC; # % atm CO2 (ppm)

# %We now need to compute the number of numerical time-steps needed to numerically solve the soil CO2 
# %partial differential equation.  The formula given from Haberman, R. (1998). Elementary applied 
# %partial differential equations (Third Edition), Prentice Hall Englewood Cliffs, NJ. (pg. 240) is
# %Ndt = (dt*Dgs)/((dz^2)*0.5).  However, we need to do max(max(Dgs)) because Dgs is a matrix.  The 
# %-0.05 # %ensures that we have more numerical time-steps than the minimum, as a precaution.
Ndt <- ceil(dt*max(max(Dgs))/((dz^2)*(0.5-0.05)));



########################################################
########################################################
########################################################

print("Running Submodels for microbial and root respiration")
# %Rall = Total microbial + root respiration at all depths and times. (mg C cm-3 hr-1):
# %Note that Rall goes from depth z = 1cm to 100cm, so we add extra column of
# %zeros (lines 91 to 93) to indicated that respiration = 0 at soil/atmosphere interface (z=0).
#function Rall =  SourceTerm(params,params_lag,Xrel,SWC_IN,SoilT_IN,SWCant_Rm_IN,SWCant_Rr_IN,SoilTant_IN,Gness_IN)

# % Submodels for microbial and root respiration.
# %
# % INPUTS:
# % params = [lagSWC_Rm lagSoilT Cstar Mstar Rmbase alpha1 alpha2 alpha3 Eo_Rm alpha4 To Km CUE lagSWC_Rr Rstar Rrbase Eo_Rr]
# % inputs = [BD PD tau solfrac Dliq]
# % Xrel  = [Rrel; Mrel], Rrel Mrel, each 1 X Nz
# % SWC_IN(t,z) = SWC (m3/m3) at time t, depth z; Nt X Nz
# % SoilT_IN(t,z) = soil temp (C) at time t, depth z; Nt X Nz
# % SWCant_IN(t,z,lagSWC) = SWC (m3/m3) at time t, depth z and lag lagSWC; Nt X Nz X Nlag 
# % SoilTant_IN(t,z,Nlag) = soil temp (C) at time t, depth z and lag lagSWC; Nt X Nz X Nlag
# % Gness_IN(t,1) = Vegetation greenness (%) at time t.  Nt X 1.
# %
# % OUTPUTS:
# % Rall.R = R (total CO2 production; mgC cm-3 hr-1).
# % Rall.Rm = Rm (microbial CO2 production; mgC cm-3 hr-1).
# % Rall.Rr = Rr (root CO2 production; mgC cm-3 hr-1)
#MATLAB#[Nt Nz]=size(SWC_IN);  #%Nt=183 and Nz=100 for daily time-step runs.
# Nt <- nrow(SWC); Nz <- ncol(SWC_IN);



#%Select the time length and depth length required from the SWC and SoilT datasets and 
#%mean-center them:
SWC <- SWC[1:Nt,1:Nz];  
SWCmean<- mean(as.matrix(SWC[1:Nt,1:Nz])); 
SWC_cent <- SWC - SWCmean;  #%mean centered SWC;
SoilT <- SoilT[1:Nt,1:Nz]+273.15; #%soil temp. in Kelvin.
SoilTmean <- mean(as.matrix(SoilT));
SoilT_cent <- SoilT - SoilTmean;  #%mean centered T;

# %Calculate the antecedent covariates by selecting the part of the dataset that corresponds 
# %to the day (for soilTant and SWCant_Rm) or week (for SWCant_Rr) in the past from datasets 
# %that we require. 
SWCant_Rm_IN <- SWCant_Rm
SWCant_Rr_IN <- SWCant_Rr
SoilTant_IN <- SoilTant

SWCant_Rm <-  sample(0,73200,replace = TRUE)
dim(SWCant_Rm) <- c(732,100)
for(L in 1:4) 
  {
    SWCant_Rm <- SWCant_Rm + (lagSWC_Rm[L] * squeeze(SWCant_Rm_IN[1:Nt,1:Nz,L]));
}


SWCant_Rr <-  sample(0,73200,replace = TRUE)
dim(SWCant_Rr) <- c(732,100)
for(L in 1:4) 
{
  SWCant_Rr <- SWCant_Rr + (lagSWC_Rr[L] * squeeze(SWCant_Rr_IN[1:Nt,1:Nz,L]));
}

SoilTant_RrRm <-  sample(0,73200,replace = TRUE)
dim(SoilTant_RrRm) <- c(732,100)
for(L in 1:4) 
  {
    SoilTant_RrRm <- SoilTant_RrRm+(lagSoilT[L]*squeeze(SoilTant_IN[1:Nt,1:Nz,L]));
  }

#%Mean center the antecedent datasets:
SWCant_Rm_cent <- SWCant_Rm-SWCmean;  
SWCant_Rr_cent <- SWCant_Rr-SWCmean;  
SoilTant_cent <- (SoilTant_RrRm+273.15)-SoilTmean;  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % MICROBIAL RESPIRATION SUBMODEL.  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % We use an extension to the DAMM model (but without O2 limitation) which incorporates 
# % microbial biomass and carbon use efficiency.  See Todd-Brown, Katherine EO, et al. "A 
# % framework for representing microbial decomposition in coupled climate models." 
# % Biogeochemistry 109.1-3 (2012): 19-33.

# % First, create matrix of tau and To values for each time and depth in model run:
# Tref=repmat(tau,Nt,Nz);    
# To=repmat(To,Nt,Nz);   


Tref <- repmat(tau,Nt,Nz);    
To <- repmat(To,Nt,Nz);   

# % Next, determine the matrix of soil organic C for each depth and time
# % point.  The vectors C and Cmic are depth dependent distributions of soil organic C and 
# % microbial C.  The formula for Cso (a Nt*Nz matrix) is taken from teh DAMM paper. We 
# % assume the distribution of microbial C is invariant with time and so copy it for each 
# % time point.  
C <- Cstar*Crel;
Cso <- repmat(as.matrix(C),Nt,1)*solfrac*(SWC^3)*Dliq;   # % always > 0
Cmic <- repmat(as.matrix(Mstar*Mrel),Nt,1);

#%Calculate Vmax.  For the first line, we need to devide the RmBase by the 
#%product of the average Cmic in teh top 10cm of soil and (1-CUE) because Vmax in
#%Todd-Brown model has different interpretation than Vmax in DAMM model.
#%Note that RmBase is the base respiration rate at Tref (since it is taken

#%from teh posterior results of the DAMM model with the PHACE data
#%(Fgrams_NEW from TG analysis).  However if we were using the alpha_s term
#%as the base rate (as done in the DAMM paper) then we would need to we
#%would need to use the formula RmBase = alpha_s/(exp(-Eas/R*Tref)).  
M=Mstar*Mrel;
RmBase_star = RmBase/mean(as.matrix(M[1:10]));  #%used for v7d onwards.

logRmB = log(RmBase_star) + (alpha1_Rm*SWC_cent) + (alpha2_Rm*SWCant_Rm_cent) + (alpha3_Rm*SWC_cent*SWCant_Rm_cent);
Eo = Eo_Rm + (alpha4*SoilTant_cent);
logVmax = logRmB + (Eo*((1/(Tref - To)) - (1/(SoilT - To))));	
Vmax = exp(logVmax);
Km=repmat(Km,Nt,Nz);

#%Bring all the diff. parts of hte model together. Computationally it would
#%be better to put the equation below on teh log-scale (and put its
#%components on log-scale as define above), but for purposes of
#%understanding, I have left it on teh raw scale.  
Rm = (Vmax*Cmic*Cso*(1-CUE))/(Km+Cso); # %units: mgC cm-3 hr-1.



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % MODEL FOR ROOT RESPIRATION PREDICTION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%Simple model of root respiration where the base respiration rate is limited by SWC and soil 
#%temp.  Units of Rr: mgC cm-3 hr-1.
Gness_IN <- Gness
fGness1=Gness_IN/max(Gness_IN);  
fGness2=1+(Gness_IN-mean(as.matrix(Gness_IN)))/(max(Gness_IN)-min(Gness_IN)); #%used for v9g
Rstar_new=Rstar*repmat(as.matrix(fGness2),1,Nz);
Cr = Rstar_new*repmat(as.matrix(Rrel),Nt,1);
logRrB = log(RrBase) + (alpha1_Rr*SWC_cent) + (alpha2_Rr*SWCant_Rr_cent) + (alpha3_Rr*SWC_cent*SWCant_Rr_cent);
Eo = Eo_Rr + (alpha4*SoilTant_cent);
logRrB_star = logRrB + (Eo*((1/(Tref-To))-(1/(SoilT-To))));
RrBase_star = exp(logRrB_star);
Rr = RrBase_star*repmat(as.matrix(fGness1),1,Nz)*Cr;  

zero_column <-  sample(0,732,replace = TRUE)
dim(zero_column) <- c(732,1)

Rr <- cbind(zero_column,Rr)
Rm <- cbind(zero_column,Rm)

# % Total soil CO2 production via respiration:
# % Total respiration (units of R are mgC cm-3 hr-1):
R <- (Rm + Rr)

########################################################
########################################################

print("Running Detect Core")

# % Convert to total source flux (mg CO2 m-3 hr-1) for each source and for each isotope 
# % (13C and 12C): Current units for R are mg C cm-3 hr-1 per layer.  So, just need to (1) 
# % convert from cm3 to m3 (mult. by 100^3); (2) convert from C to CO2 (mult. by 44/12).

#print("Calculating S")
S <- R*(44/12)*(100^3);

#print("Calculating Sr")
Sr <- Rr*(44/12)*(100^3);

#print("Calculating Sm")
Sm <- Rm*(44/12)*(100^3);

# %Partition the soil CO2 sources into the isotopic components:
# %The lines below (115-120) are based on solving eqns (1) and (2) simultaneously for 
# %[13C] (=S13root) and [12C] (=S12root):
# %delta13Cr = [(Rr – Rstd)/Rstd] * 1000   ------------------- (1)
# %where Rr = [13C] / [12C]
# %[13C] + [12C] = Sr    ------------------------------------- (2)
# %We solve for S13micr and S12micr in the same way.
# %Note that the d13Cr and d13Cm refer to the isotopic signature of the
# %micr. and root respiration at each depth, i.e. it's not hte isotopic
# %signature of the surface respired CO2.  For the 2nd DETECT paper, we
# %may wish to allow these d13Cr and d13Cm terms to vary with time (but 
# %still keep fixed with depth).  
Rstd <- 0.0112372;  # % Rstd for 13C/12C (Pee Dee Belemite):

d13Cr_fixed <- (as.numeric(unlist(d13Cr, use.names=FALSE)))
denomR <- (1000+(1000+repmat(d13Cr_fixed,Nt+1))*Rstd);

#print("Calculating d13Cr")
d13Cm_fixed <- (as.numeric(unlist(d13Cm, use.names=FALSE)))
#print("Calculating denomM")
denomM = (1000+(1000+repmat(d13Cm_fixed,Nt+1))*Rstd);

#print("Calculating Components of S type")
S13root = Sr*(1000+repmat(d13Cr_fixed,Nt+1))*Rstd/denomR;
S12root = 1000*Sr/denomR;
S13micr = Sm*(1000+repmat(d13Cm_fixed,Nt+1))*Rstd/denomM;
S12micr = 1000*Sm/denomM;

#print("Creating Stype Array")
Stype <- as.array(array(rep(NaN, 4*732*101), dim = c(4, 732, 101)))  
for (i in 1:732)
  {
    for (j in 1:100)
    {
      Stype[1,i,j] <- S13root[i,j]
      Stype[2,i,j] <- S12root[i,j]
      Stype[3,i,j] <- S13micr[i,j]
      Stype[4,i,j] <- S12micr[i,j]
      
    }
  }

# %Convert upper boundary condition (upbc) (ie atmospheric [CO2] or [CO2] at the soil
# %surface) from ppm units to mg CO2 m-3 units.  Note that 101300 is the standard 
# %pressure in Pa, R=8.134 is the gas constant, and the soil temp T is in Kelvin.
#print("Calculating upbc")
upbc = (catm*44/1000)*101300./(8.314*T_T[,1]); 

# %The Partial differential equation (equation 1 of Ryan et al., in review) is solved 
# %four times for different types of soil CO2: (1) the d12C of soil CO2 from microbial 
# %sources; (2) the d13C of soil CO2 from microbial sources; (3) the d12C of soil CO2 
# %from root sources; (4) the d13C of soil CO2 from root sources.  This for
# %loop computes the initial and boundary conditions for the four types of
# %soil CO2.
#print("Creating cctype matrix")
cctype <- as.array(array(rep(-9999, 4*732*101), dim = c(4, Nz+1, Nt+1)))  
fluxtype <- zeros(4, Nt);
for (i in 1:4)
  {
    cctype[i,2:Nz,1] <- as.array(as.matrix((squeeze(Stype[i,1,2:Nz])/S[1,2:Nz])*as.numeric(ic[2:Nz])));
    cctype[i,Nz+1,1] <- cctype[i,Nz,1];
    cctype[i,1,1] <- as.array(as.matrix((squeeze(Stype[i,1,2])/S[1,2])*as.numeric(ic[1])));
    cctype[i,1,2:(Nt+1)] <- as.array(as.matrix((squeeze(Stype[i,1,2])/S[1,2])*upbc));
  }

# %Run main For Loop.  This computes the soil CO2 concentration for each of the four 
# %types described above by numerically solving the partial differential equation.
ctype <- zeros(4,Nz+1);
dt_old <- dt;
dt <- dt_old/Ndt;
Nt_run <- as.numeric(Nt_run)


print("Running main DETECT Loop")
Source_type <- c("S13root","S12root","S13micr","S12micr")
Start_time <- Sys.time()
for (t in 1:Nt_run) {
  New_DGS_1_3 <- as.numeric(c(Dgs_13C[t,1],Dgs_13C[t, ]))
  New_DGS_2_4 <- as.numeric(c(Dgs_12C[t,1],Dgs_12C[t, ]))
  print("Current Iteration:")
  print(t)
  if (t > 1) {
    for (i in 1:4) {
      cctype[i,1,t] <- as.array(as.matrix((squeeze(cctype[i,2,t-1])/sum(cctype[ ,2,t-1]))*upbc[t-1,1]));
      cctype[i,Nz+1,t] <- cctype[i,Nz,t];
    }
  }
  ctype = squeeze(cctype[ , ,t]);
  for (tt in 1:Ndt){
    
    ctype_old <- ctype;
    for (z in 2:Nz) {
      for (i in 1:4) {
        if (i == 1 || i == 3) {
          Dgs_t <- New_DGS_1_3
        } 
        else {
          Dgs_t <- New_DGS_2_4
        }
        
        if (z==Nz) {
          ctype[i,z] = ctype_old[i,z] + dt*((Dgs_t[z]/(dz*dz))*(ctype_old[i,z-1]-ctype_old[i,z])+ Stype[i,t,z]);
        }
        else {
          ctype[i,z] = ctype_old[i,z] + dt*((Dgs_t[z]/(dz*dz))*(ctype_old[i,z+1]-2*ctype_old[i,z]+ctype_old[i,z-1]) + ((Dgs_t[z+1]-Dgs_t[z-1])*(ctype_old[i,z+1]-ctype_old[i,z-1]))/(4*dz^2) + Stype[i,t,z]);
        }
      }
    }
  }
  cctype[ , ,t+1] <- c(ctype[1:4,1:Nz], ctype[1:4,Nz]);
  for (i in 1:4){
    if ((i==1) || (i==3)){
      MeanDgs <- Dgs_13C[t,1];
    }
    else {
      MeanDgs <- Dgs_12C[t,1];
    }
    fluxtype[i,t] = (MeanDgs*(cctype[i,2,t]-cctype[i,1,t])/dz)/(44*3.6);
  }
  # %dividing by (44*3.6) converts units of soil CO2 from mg CO2 m-2 hr-1 to
  # %umol CO2 m-2 s-1.
}

End_time <- Sys.time()
Time_to_run <- End_time - Start_time; print(Time_to_run)

writeMat(con= paste0("./Outputs/", format(End_time, "%d_%b_%Y_%H_%M"), "_cctype.mat"), x=(cctype))
# writeMat(con= paste0("./output/", format(End_time, "%d_%b_%Y_%H_%M"), "_cctype.mat"), x=(cctype))
# writeMat(con="./output/Dgs.mat", x=(Dgs))
# writeMat(con="./output/phig.mat", x=(phig))
# 
# writeMat(con="./output/R.mat", x=(R))
# writeMat(con="./output/S.mat", x=(S))
# writeMat(con="./output/Rm.mat", x=(Rm))
# writeMat(con="./output/Rr.mat", x=(Rr))
# writeMat(con="./output/Stype.mat", x=(Stype))
# writeMat(con="./output/Dgs_t.mat", x=(Dgs_t))
