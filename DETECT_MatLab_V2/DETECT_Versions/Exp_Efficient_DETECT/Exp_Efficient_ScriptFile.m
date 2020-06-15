clear all
close all

% Run model with plants roots removed from soil.
Bare_Soil = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SCRIPT FILE FOR RUNNING THE SOIL CO2 TRANSPORT AND PRODUCTION MODEL (a.k.a. the DETECT MODEL)
% ---------------------------------------------------------------------------------------------
%
% This model is described in the paper: 
% Ryan, E. M., K. Ogle, H. Kropp, Y. Carrillo, K. E. Samuels-Crow, E. Pendall (in review) Modelling the
% soil CO2 transport and production with non-scalar diffusion and source terms, part 1: model description 
% and testing the steady-steady assumption. Submitted to the Journal Of Advances in Modelling Earth 
% Systems. 
%
% The purpose of this script file is threefold: (1) to read in the 'Data_6hourly.mat' file where all the 
% data needed to drive the model is stored; (2) to define the values of the model parameters; (3) to run 
% the DETECT model which is coded up in a separate Matlab file called DETECT_July2017.m.  The other Matlab 
% file ('SourceTerm_July2017.m') contains the code for calculating the production term (S term in Ryan et 
% al. paper) and is read in by DETECT_July2017.m.
%
% If you are using this for the first time and wish to try it out, please press the 'Run' buttom on the 
% top panel.  You can plot a time-series of different model outputs (e.g. soil respiration) using the 
% Plot_output.m file.
%
% If you wish to run this model based on driving data from your site, please see the manual to make sure 
% you have all the site data required.  Please also see the 'Data_6hourly_Creator.m' file, which creates 
% the .mat file needed in order to run this script file. Finally, you will also need estimates of the 
% model parameters which are described below.  
%
% This version of the script file is the same as the working version used in the Ryan et al. paper except 
% that (1) the code has been tidied up and (2) BD, b and SWC100 (parameters for the soil water potential 
% function) are now depth dependent.                                                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Things to set before running this script file (for details, see Ryan et al, in review):
running_antecedent_vers=0;  % Including (=1) or not (=0) the antecedent variables in the DETECT model.
runmodel=1;                 % Whether to run the model (=1) or not(=0).
Nt_run=183*4;               % The number of time steps the user wishes to run the DETECT model for. For 
                            % Ryan et al. (in review), max(Nt_run)=183*4, where 183 is the no. of days 
                            % model is run for and 4 is the no. of 6 hourly time-steps each day.
Nt_beg=92;                  % First DoY to use for computer experiment.
Nt_end=274;                 % Final DoY to use for computer experiment.
writefile=1;                % Whether to write the output to mat and csv files (=1) or not (=0).

% Read in the driving data and site data.  For Ryan et al. (in review), SWC (Soil Water Content) and SoilT 
% (Soil Temperature) are output from another model called HYDRUS which estimates these variables based on 
% measurements of SWC and SoilT at three depths at the site.  The AtmPress and Gness data are measured at 
% site.  Temporal gap-filling of these four sets of measurements were made
% for times of missing data.  The Nt and Nz refer to the number of time 
load(['Exp_Efficient_Inputs/PHACEdata_6hourly.mat'])  %Control treatment plot from PHACE site for 2008
SWC = Data.SWC;                %Nt * Nz
SoilT = Data.SoilT;            %Nt * Nz
SWCant_Rm = Data.SWCantRm;     %Nt * Nz * Nlag  (Nlag = 4 weeks in the past)
SWCant_Rr = Data.SWCantRr;     %Nt * Nz * Nlag  (Nlag = 4 weeks in the past)
SoilTant = Data.SoilTant;      %Nt * Nz * Nlag  (Nlag = 7 days in the past)
Gness = Data.Gness;            %Nt * 1
AtmPress = Data.AtmPress;      %Nt * 1
params_swp = Data.params_swp;  %Nz * 3
Xrel = Data.Xrel;              %Nz * 3
ic = Data.cinit;               %1 * Nz+1
d13Cr = Data.d13Cr;            %1 * Nz+1
d13Cm = Data.d13Cm;            %1 * Nz
dx = Data.dx;                  %1 * 2
BC = Data.BC;                  %Nt * 1    


% Parameters used by the DETECT model.  For a detailed description of the parameters, please see Ryan et al. 
% (in review).  Parameter values used are mostly taken from the following papers:
% Ryan et al. (2015) Antecedent moisture and temperature conditions modulate the response of ecosystem 
% respiration to elevated CO2 and warming
% Davidson et al. (2012) The Dual Arrhenius and Michaelis–Menten kinetics model for decomposition of soil 
% organic matter at hourly to seasonal time scales Tucker et al. (2012) Does declining carbon-use efficiency 
% explain thermal acclimation of soil respiration with warming?

%Additional parameters required for the Rr model:
params(1) = 111.54;  % Rstar = total root biomass (mgC) mass in a one cm2 column of soil 1m deep. 
if Bare_Soil == 1,
    params(1) = 0;
end
params(2) = 0.00006; % RrBase, mass-specific root respiration base rate at a tau degC (mgC cm-3 hr-1).
params(3) = 11.65;    % alpha1 for Rr; SWC param.    
params(4) = 20.72;   % alpha2 for Rr; The main effect of SWCant on the log-scaled Rrbase
params(5) = -164.16; % alpha3 for Rr; The main effect of SWC*SWCant on the log-scaled Rrbase

%Parameters used by the Rm submodel.
params(6) = 711.6094; % Cstar = total soil C (mgC) mass in a one cm2 column of soil 1m deep. 
params(7) = 12.2568;  % Mstar = total MBC (mgC) mass in a one cm2 column of soil 1m deep. 
params(8) = 0.0015;   % RmBase; (mg C cm-3 hr-1); based on value used in Davidson et al. (2012)
params(9) = 14.05;    % alpha1 (mg C cm-3 hr-1); based on value used in Ryan et al. (2015) 
params(10) = 11.05;    % alpha2; The main effect of SWCant on the log-scaled Rmbase;  
params(11) = -87.55;   % alpha3; The main effect of SWC*SWCant on the log-scaled Rmbase; 
params(12) = 0.00001; % Km
params(13) = 0.8;     % CUE; Carbon Use Efficiency (mgC/mgC).  Taken from Tucker et al. (2012) 
params(14) = 0.0004;  % solfrac; Fraction of substrate that Sx that is soluble (called 'p' in DAMM paper). 
params(15) = 3.17;    % Dliq; Diffusivity of substrate Sx that is liquid.  Dimensionless.
                 
params(16) = 324.6;    % Eo_Rm; Eo term used in Rm submodel (Kelvin). Use value from Ryan et al. (2015)
params(17) = 324.6;   % Eo_Rr; Eo term used in Rr submodel (Kelvin); Use value from Ryan et al. (2015)
params(18) = 227.5;   % To; Temp. parameter used in the Lloyd & Taylor temp. (Kelvin)
params(19) = -4.67;    % alpha4; The main effect of SoilTant on the log-scaled Rmbase;
params(20) = 10;      % tau; Base temperature in Kelvin.  We use 10 degC.
                                          
if (running_antecedent_vers==0)
    params(4) = 0; params(5) = 0; params(10) = 0; params(11) = 0; params(19) = 0; 
    ant=[''];
else
    ant=['_ant'];
end

%Lag parameters.  The first row describes the importance of each of hte past
%4 days on Rm in Source term submodel.  The 2nd row describes the importance 
%of each of hte past 4 weeks on Rr in Source term submodel.
params_lag = [0.75 0.25 0 0; 
              0.2 0.6 0.2 0;
              0.25 0.25 0.25 0.25];  

% Run the soil CO2 flux and concentration model.
if (runmodel==1)
    t1=clock;    %start time for model run.
    [out info] = Exp_Efficient_DETECTmodel(Nt_run,SWC,SoilT,SWCant_Rm,SWCant_Rr,SoilTant,Gness,AtmPress,params_swp,Xrel,ic,d13Cr,d13Cm,dx,BC,params,params_lag);
    t2=clock;
end

%Out is a structure that contains the soil [CO2] at different depths for
%each time point (out.CO2); the soil [CO2] contributed by different
%sources by isotope (12C and 13C) (out.CO2type); soil CO2 flux at the surface
%(out.CO2flux); soil CO2 flux at the surface partitioned into
%contributions by microbes and roots by isotope (out.CO2fluxtype); and the
%13C signature soil respiration (out.d13C).

%Save the output to mat files. 
%We also save the out.CO2 matrix to a csv file to view it easier.
DayN=[(Nt_beg-1):0.25:Nt_end]';
if (writefile==1)
    save([pwd '/DETECT_Versions/Exp_Efficient_DETECT/Exp_Efficient_Outputs/output_6hourly' num2str(ant) '.mat'],'out')
    save([pwd '/DETECT_Versions/Exp_Efficient_DETECT/Exp_Efficient_Outputs/info_6hourly' num2str(ant) '.mat'],'info')
    myfile4 = fopen([pwd '/DETECT_Versions/Exp_Efficient_DETECT/Exp_Efficient_Outputs/outCO2_6hourly' num2str(ant) '.csv'] ,'w+');
    fprintf(myfile4,['DayN,' 'd0,' 'd1,' 'd2,' 'd3,' 'd4,' 'd5,' 'd6,' 'd7,' 'd8,' 'd9,' 'd10,' 'd11,' 'd12,' 'd13,' 'd14,' 'd15,' 'd16,' 'd17,' 'd18,' 'd19,' 'd20,' 'd21,' 'd22,' 'd23,' 'd24,' 'd25,' 'd26,' 'd27,' 'd28,' 'd29,' 'd30,' 'd31,' 'd32,' 'd33,' 'd34,' 'd35,' 'd36,' 'd37,' 'd38,' 'd39,' 'd40,'  'd41,' 'd42,' 'd43,' 'd44,' 'd45,' 'd46,' 'd47,' 'd48,' 'd49,' 'd50,'  'd51,' 'd52,' 'd53,' 'd54,' 'd55,' 'd56,' 'd57,' 'd58,' 'd59,' 'd60,'  'd61,' 'd62,' 'd63,' 'd64,' 'd65,' 'd66,' 'd67,' 'd68,' 'd69,' 'd70,'  'd71,' 'd72,' 'd73,' 'd74,' 'd75,' 'd76,' 'd77,' 'd78,' 'd79,' 'd80,'  'd81,' 'd82,' 'd83,' 'd84,' 'd85,' 'd86,' 'd87,' 'd88,' 'd89,' 'd90,'  'd91,' 'd92,' 'd93,' 'd94,' 'd95,' 'd96,' 'd97,' 'd98,' 'd99,' 'd100,\n']);
    fclose(myfile4);
    dlmwrite([pwd '/DETECT_Versions/Exp_Efficient_DETECT/Exp_Efficient_Outputs/outCO2_6hourly' num2str(ant) '.csv'],[DayN out.CO2], '-append','precision', '%.3f');
end
