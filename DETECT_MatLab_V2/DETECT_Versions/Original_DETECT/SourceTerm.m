function Rall =  SourceTerm(params,params_lag,Xrel,SWC_IN,SoilT_IN,SWCant_Rm_IN,SWCant_Rr_IN,SoilTant_IN,Gness_IN)
% Submodels for microbial and root respiration.
%
% INPUTS:
% params = [lagSWC_Rm lagSoilT Cstar Mstar Rmbase alpha1 alpha2 alpha3 Eo_Rm alpha4 To Km CUE lagSWC_Rr Rstar Rrbase Eo_Rr]
% inputs = [BD PD tau solfrac Dliq]
% Xrel  = [Rrel; Mrel], Rrel Mrel, each 1 X Nz
% SWC_IN(t,z) = SWC (m3/m3) at time t, depth z; Nt X Nz
% SoilT_IN(t,z) = soil temp (C) at time t, depth z; Nt X Nz
% SWCant_IN(t,z,lagSWC) = SWC (m3/m3) at time t, depth z and lag lagSWC; Nt X Nz X Nlag 
% SoilTant_IN(t,z,Nlag) = soil temp (C) at time t, depth z and lag lagSWC; Nt X Nz X Nlag
% Gness_IN(t,1) = Vegetation greenness (%) at time t.  Nt X 1.
%
% OUTPUTS:
% Rall.R = R (total CO2 production; mgC cm-3 hr-1).
% Rall.Rm = Rm (microbial CO2 production; mgC cm-3 hr-1).
% Rall.Rr = Rr (root CO2 production; mgC cm-3 hr-1)

[Nt Nz]=size(SWC_IN);  %Nt=183 and Nz=100 for daily time-step runs.

%Parameters used by the Rm model:
Rstar = params(1);     % Rstar = sum(Rrel) = total root biomass (mgC) mass in a one cm2 column of soil 1m deep. 
RrBase = params(2);    % RrBase, mass-specific root respiration base rate at a tau degC (mgC cm-3 hr-1).
alpha1_Rr = params(3); % The main effect of SWC on the log-scaled Rrbase; 
alpha2_Rr = params(4); % The main effect of SWCant on the log-scaled Rrbase; 
alpha3_Rr = params(5); % The main effect of SWC*SWCant on the log-scaled Rrbase; 

Cstar = params(6);     % Cstar in M-M model = total soil C (mgC) mass in a one cm2 column of soil 1m deep. 
Mstar = params(7);     % Mstar in M-M model = total MBC (mgC) mass in a one cm2 column of soil 1m deep. 
RmBase = params(8);    % base respiration rate for microbial respiration function at tau deg. C (mgC cm-3 hr-1).
alpha1_Rm = params(9); % The main effect of SWC on the log-scaled Rmbase; 
alpha2_Rm = params(10); % The main effect of SWCant on the log-scaled Rmbase; 
alpha3_Rm = params(11); % The main effect of SWC*SWCant on the log-scaled Rmbase; 
Km = params(12);       % Michaelis constants (?) from DAMM paper
CUE = params(13);      % Carbon Use Efficiency (mgC/mgC)
solfrac = params(14);    % Fract. of substrate that Sx that is soluble (called 'p' in DAMM paper). Dimesionless.
Dliq = params(15);       % Diffusivity of substrate Sx that is liquid.  Dimensionless.

Eo_Rm = params(16);     % Eo term used in Rm submodel (Kelvin)  
Eo_Rr = params(17);     % Eo term used in Rr submodel (Kelvin)  
To = params(18);       % Temp. parameter used in the Lloyd & Taylor temp. response type function for Rm (Kelvin)
alpha4 = params(19);    % The main effect of SoilTant on the log-scaled Rmbase Rrbase; 
tau = params(20)+273.15; % Base temperature in Kelvin.

%Lag parameters
lagSWC_Rm=params_lag(1,:);
lagSWC_Rr=params_lag(2,:);
lagSoilT=params_lag(3,:);

%The vectors of the relative depth-dependent distributions of microbial biomass
%carbon (Mrel), total root biomass (Rrel) and soil organic carbon (Crel).
Crel = Xrel(1,:);
Mrel = Xrel(2,:);
Rrel = Xrel(3,:);

%Select the time length and depth length required from the SWC and SoilT datasets and 
%mean-center them:
SWC=SWC_IN(1:Nt,1:Nz);     
SWCmean=mean(mean(SWC));
SWC_cent=SWC-SWCmean;  %mean centered SWC;
SoilT=SoilT_IN(1:Nt,1:Nz)+273.15; %soil temp. in Kelvin.
SoilTmean=mean(mean(SoilT));
SoilT_cent=SoilT-SoilTmean;  %mean centered T;

%Calculate the antecedent covariates by selecting the part of the dataset that corresponds 
%to the day (for soilTant and SWCant_Rm) or week (for SWCant_Rr) in the past from datasets 
%that we require. 
SWCant_Rm=zeros(Nt,Nz);
for L=1:4,
    SWCant_Rm=SWCant_Rm+(lagSWC_Rm(1,L)*squeeze(SWCant_Rm_IN(1:Nt,1:Nz,L)));
end

SWCant_Rr=zeros(Nt,Nz);
for L=1:4,
    SWCant_Rr=SWCant_Rr+(lagSWC_Rr(1,L)*squeeze(SWCant_Rr_IN(1:Nt,1:Nz,L)));
end

SoilTant_RrRm=zeros(Nt,Nz);
for L=1:4,
    SoilTant_RrRm=SoilTant_RrRm+(lagSoilT(1,L)*squeeze(SoilTant_IN(1:Nt,1:Nz,L)));
end

%Mean center the antecedent datasets:
SWCant_Rm_cent=SWCant_Rm-SWCmean;  
SWCant_Rr_cent=SWCant_Rr-SWCmean;  
SoilTant_cent=(SoilTant_RrRm+273.15)-SoilTmean;  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MICROBIAL RESPIRATION SUBMODEL.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We use an extension to the DAMM model (but without O2 limitation) which incorporates 
% microbial biomass and carbon use efficiency.  See Todd-Brown, Katherine EO, et al. "A 
% framework for representing microbial decomposition in coupled climate models." 
% Biogeochemistry 109.1-3 (2012): 19-33.

% First, create matrix of tau and To values for each time and depth in model run:
Tref=repmat(tau,Nt,Nz);    
To=repmat(To,Nt,Nz);   

% Next, determine the matrix of soil organic C for each depth and time
% point.  The vectors C and Cmic are depth dependent distributions of soil organic C and 
% microbial C.  The formula for Cso (a Nt*Nz matrix) is taken from teh DAMM paper. We 
% assume the distribution of microbial C is invariant with time and so copy it for each 
% time point.  
C = Cstar*Crel;
Cso = repmat(C,Nt,1).*solfrac.*(SWC.^3).*Dliq;  % always > 0
Cmic=repmat(Mstar*Mrel,Nt,1);

%Calculate Vmax.  For the first line, we need to devide the RmBase by the 
%product of the average Cmic in teh top 10cm of soil and (1-CUE) because Vmax in
%Todd-Brown model has different interpretation than Vmax in DAMM model.
%Note that RmBase is the base respiration rate at Tref (since it is taken
%from teh posterior results of the DAMM model with the PHACE data
%(Fgrams_NEW from TG analysis).  However if we were using the alpha_s term
%as the base rate (as done in the DAMM paper) then we would need to we
%would need to use the formula RmBase = alpha_s/(exp(-Eas/R*Tref)).  
M=Mstar*Mrel;
RmBase_star = RmBase/mean(M(1:10));  %used for v7d onwards.
logRmB = log(RmBase_star) + (alpha1_Rm*SWC_cent) + (alpha2_Rm*SWCant_Rm_cent) + (alpha3_Rm*SWC_cent.*SWCant_Rm_cent);
Eo = Eo_Rm + (alpha4*SoilTant_cent);
logVmax = logRmB + (Eo.*((1./(Tref - To)) - (1./(SoilT - To))));	
Vmax = exp(logVmax);
Km=repmat(Km,Nt,Nz);

%Bring all the diff. parts of hte model together. Computationally it would
%be better to put the equation below on teh log-scale (and put its
%components on log-scale as define above), but for purposes of
%understanding, I have left it on teh raw scale.  
Rm = (Vmax.*Cmic.*Cso.*(1-CUE))./(Km+Cso);  %units: mgC cm-3 hr-1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL FOR ROOT RESPIRATION PREDICTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simple model of root respiration where the base respiration rate is limited by SWC and soil 
%temp.  Units of Rr: mgC cm-3 hr-1.
fGness1=Gness_IN/max(Gness_IN);  
fGness2=1+(Gness_IN-mean(Gness_IN))./(max(Gness_IN)-min(Gness_IN)); %used for v9g
Rstar_new=Rstar*repmat(fGness2,1,Nz);
Cr = Rstar_new.*repmat(Rrel,Nt,1);
logRrB = log(RrBase) + (alpha1_Rr*SWC_cent) + (alpha2_Rr*SWCant_Rr_cent) + (alpha3_Rr*SWC_cent.*SWCant_Rr_cent);
Eo = Eo_Rr + (alpha4*SoilTant_cent);
logRrB_star = logRrB + (Eo.*((1./(Tref-To))-(1./(SoilT-To))));
RrBase_star = exp(logRrB_star);
Rr = RrBase_star.*repmat(fGness1,1,Nz).*Cr;  

% Total soil CO2 production via respiration:
% Total respiration (units of R are mgC cm-3 hr-1):
R = (Rm + Rr);

Rall.R = R; Rall.Rm= Rm; Rall.Rr= Rr;