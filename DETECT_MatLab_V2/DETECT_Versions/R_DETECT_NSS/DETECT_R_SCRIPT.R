print("Initializing script")
print("###################")
print("Clearing enviroment")
rm(list=ls());print("Enviroment is cleared");
print("###################")
print("Adding required libraries")

packages <- c("rmatio", "R.matlab", "pracma", "tidyverse", "imager", "here")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
library(rmatio)
library(R.matlab)
library(pracma)
library(tidyverse)
library(imager)

print("###################")
print("Setting up working directory")
setwd(paste(here:::here(),'/DETECT_Versions/R_DETECT_NSS/',sep = ''));print("Working directory is changed to:"); getwd()
print("###################")
# 
# print("###################")
# print("###################")

source("ReadInData_1.R")

# %# %# %# %# %# %# %# %# %# %# %# %#
print("Core code")
# %# %# %# %# %# %# %# %# %# %# %# %# 

# % Based on soil CO2 production model developed by Fang & Moncrieff (1999)
# % Agricultural and Forest Meteorology 95:225-236. But, simplify model to
# % only deal with gas phase CO2 as is done in Hui & Luo (2004) Global
# % Biogeochemical Cycles Vol 18:GB4029.
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

# %Rall = Total microbial + root respiration at all depths and times. (mg C cm-3 hr-1):
# %Note that Rall goes from depth z = 1cm to 100cm, so we add extra column of
# %zeros (lines 91 to 93) to indicated that respiration = 0 at soil/atmosphere interface (z=0).
source("ReadInSourceTerm_SOIL_MICRO.R")

# % Convert to total source flux (mg CO2 m-3 hr-1) for each source and for each isotope 
# % (13C and 12C): Current units for R are mg C cm-3 hr-1 per layer.  So, just need to (1) 
# % convert from cm3 to m3 (mult. by 100^3); (2) convert from C to CO2 (mult. by 44/12).
# print("Calculating S")
S <- R*(44/12)*(100^3);
# print("Calculating Sr")
Sr <- Rr*(44/12)*(100^3);
# print("Calculating Sm")
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
# print("Calculating d13Cr")
d13Cm_fixed <- (as.numeric(unlist(d13Cm, use.names=FALSE)))
# print("Calculating denomM")
denomM = (1000+(1000+repmat(d13Cm_fixed,Nt+1))*Rstd);
# print("Calculating Components of S type")
S13root = Sr*(1000+repmat(d13Cr_fixed,Nt+1))*Rstd/denomR;
S12root = 1000*Sr/denomR;
S13micr = Sm*(1000+repmat(d13Cm_fixed,Nt+1))*Rstd/denomM;
S12micr = 1000*Sm/denomM;

Stype <- as.array(array(rep(NaN, 4*732*101), dim = c(4, 732, 101)))  
Stype[1,,] <- data.matrix(S13root, rownames.force = NA)
Stype[2,,] <- data.matrix(S12root, rownames.force = NA)
Stype[3,,] <- data.matrix(S13micr, rownames.force = NA)
Stype[4,,] <- data.matrix(S12micr, rownames.force = NA)

# %Convert upper boundary condition (upbc) (ie atmospheric [CO2] or [CO2] at the soil
# %surface) from ppm units to mg CO2 m-3 units.  Note that 101300 is the standard 
# %pressure in Pa, R=8.134 is the gas constant, and the soil temp T is in Kelvin.
# print("Calculating upbc")
upbc = (catm*44/1000)*101300./(8.314*T_T[,1]); 

# %The Partial differential equation (equation 1 of Ryan et al., in review) is solved 
# %four times for different types of soil CO2: (1) the d12C of soil CO2 from microbial 
# %sources; (2) the d13C of soil CO2 from microbial sources; (3) the d12C of soil CO2 
# %from root sources; (4) the d13C of soil CO2 from root sources.  This for
# %loop computes the initial and boundary conditions for the four types of
# %soil CO2.
# print("Creating cctype matrix")
cctype <- as.array(array(rep(-9999, 4*732*101), dim = c(4, Nz+1, Nt+1)))  
fluxtype <- zeros(4, Nt);


for (i in 1:4){
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
Source_type <- c("S13root","S12root","S13micr","S12micr")

# source("Main__Loop_Development_NSS.R")
# Source_type <- c("S13root","S12root","S13micr","S12micr")
# Start_time <- Sys.time()
for (t in 1:Nt_run) {
  
  New_DGS_1_3 <- as.numeric(c(Dgs_13C[t,1],Dgs_13C[t, ]))
  New_DGS_2_4 <- as.numeric(c(Dgs_12C[t,1],Dgs_12C[t, ]))
  # print("Current Iteration:")
  # print(t)
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
}
# End_time <- Sys.time()
# Time_to_run <- End_time - Start_time; print(Time_to_run)


writeMat(con="./Outputs/output_6hourly.mat", out.CO2type=(cctype))
# writeMat(con="./Outputs/output_6hourly.mat", Dgs=(Dgs))
# writeMat(con="./Outputs/output_6hourly.mat", phig=(phig))
# writeMat(con="./Outputs/output_6hourly.mat", R=(R))
# writeMat(con="./Outputs/output_6hourly.mat", S=(S))
# writeMat(con="./Outputs/output_6hourly.mat", Rm=(Rm))
# writeMat(con="./Outputs/output_6hourly.mat", Rr=(Rr))
# writeMat(con="./Outputs/output_6hourly.mat", Stype=(Stype))
# writeMat(con="./Outputs/output_6hourly.mat", Dgs_t=(Dgs_t))


# writeMat(con="./output/cctype.mat", CO2flux=(cctype))
# writeMat(con="./output/Dgs.mat", Dgs=(Dgs))
# writeMat(con="./output/phig.mat", phig=(phig))
# writeMat(con="./output/R.mat", R=(R))
# writeMat(con="./output/S.mat", S=(S))
# writeMat(con="./output/Rm.mat", Rm=(Rm))
# writeMat(con="./output/Rr.mat", Rr=(Rr))
# writeMat(con="./output/Stype.mat", Stype=(Stype))
# writeMat(con="./output/Dgs_t.mat", Dgs_t=(Dgs_t))
