clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This m file creates the mat file which stores all the environmental and site data required to drive the 
% DETECT model.  The information needed
% before running this m file are as follows (1-5 are read in as csv files):
% (1) soil water content data (m2/m2) for all (100) depths and days. Edit line 51.
% (2) soil temperature data (degree C) for all depths and times. Edit line 52.
% (3) air pressure data for the site for all times.  Edit line 53.
% (4) greenness data for the site for all times.  Edit line 54.
% (5) the initial soil CO2 values for all depths.  Edit line 55.
% (6) Parameters required to estimate the Soil water potential: bulk % density, b and SWC100 (Soil water content
%     at -100cmH20).  Edit line 56.
% (7) The atmospheric CO2 concentration for the site (edit line 57).
% (8) Parameter values for the gamma and exponential functions used to describe the distribution by depth of
%     microbial biomass, root biomass and soil organic matter (edit lines 132-134).  These values can be obtained 
%     by comparing different distributions (by trial and improvement) with measurements of root biomass, microbial 
%     biomass and soil organic matter from a handful of depths.  Set write.file=0 (line 31) and Xrel_params=1.
%     (line 32) when doing this comparison.  Once finished, set write.file=1 and Xrel_params=0 (edit lines 183-236).
% (9) The 13C signature of the CO2 respired by roots and microbes (edit lines 159-160).
% (10) The time and depth increments (edit lines 34-37)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Things to specify:
site=['PHACE'];              %Name of field site.
write_file = 1;              %whether to create hte mat file (=1) or not (=0).
Xrel_params = 0;             %See lines 18-22 above. 
year=2008;                   %Year chosen to run the DETECT model for.
Nt_d=183;                    %No. of days to run DETECT for.  We assume length of time-step is 6 hours.
Nz=100;                      %No. of depths.  We assume the depth increment is 1cm.
Nt_beg=92;                   %First DoY to use for experiment.
Nt_beg_ant=1;                %DoY values to start at for SWC and soilT data - we start earlier for antecedent calculations.

%Creating the time vectors.
Nt_end = Nt_beg + Nt_d - 1;
Nt=Nt_d*4;
time1_d=[Nt_beg:Nt_end]';      %DoY numbers to use for DETECT model run.
time1_6h=[(Nt_beg-0.75):0.25:Nt_end]'*4;  %d=daily; 6h=6hourly
time2_d=[1:Nt_end]';              %DoY values used for SWC and soilT data - we start earlier for antecedent calculations.
time2_6h=[0.25:0.25:Nt_end]'*4;

%Load data.  SWC and SoilT are HYDRUS output based on gapfilled data at 3
%depths.  The AtmPress and Gness data are gapfilled by NARR and Ed Ryan
%respectively. cinit is the vector of t=0 values of the soil CO2.  it is obtained 
%by running the model for 2007 and taking the soil CO2 values for 31st Dec.
SoilT_6h_IN=csvread([num2str(site) '_SoilT_2008_6hourly.csv'],1,2);
SWC_d_IN=csvread([num2str(site) '_SWC_2008_daily.csv'],1,1);
AtmPress_d_IN=csvread([num2str(site) '_Press_2008_daily.csv'],1,1);
Gness_d_IN=csvread([num2str(site) '_Gness_2008_daily.csv'],1,1);
cinit=csvread([num2str(site) '_soilCO2_initial.csv'],1,0);
params_swp=csvread([num2str(site) '_params_SWP.csv'],1,1);
atmCO2_IN=csvread([num2str(site) '_atmCO2_2008_6hourly.csv'],1,1);

%This is to make sure we haven't read in any additional rows or cols:
SoilT_6h=SoilT_6h_IN(1:1464,1:101);
SWC_d=SWC_d_IN(1:366,1:101);
AtmPress_d=AtmPress_d_IN(1:366,1);
Gness_d=Gness_d_IN(1:366,1);

%Expand the SWC, Press and Gness datasets to be 6 hourly:
SWC_6h=ones(length(SWC_d)*4,Nz+1)*-9999;
AtmPress_6h=ones(length(AtmPress_d)*4,1)*-9999;
Gness_6h=ones(length(Gness_d)*4,1)*-9999;
for i=1:size(SWC_d,1)
    endT=i*4; begT=endT-3;
    SWC_6h(begT:endT,1:101)=repmat(SWC_d(i,1:(Nz+1)),4,1);
    AtmPress_6h(begT:endT,1)=repmat(AtmPress_d(i,1),4,1);
    Gness_6h(begT:endT,1)=repmat(Gness_d(i,1),4,1);
end

%Select the days we want for the current data:
z=2:1:(Nz+1);   %Row number in matrix of the different datasets, i.e. row=2 is depth=1cm.
SoilT=SoilT_6h(time1_6h,z);
SWC=SWC_6h(time1_6h,z);
AtmPress=AtmPress_6h(time1_6h,1)/10;  %divide by 10 to convert from hPa to kPa.;
Gness=Gness_6h(time1_6h,1)/100;  %divide by 100 to convert from % to 0-1 scale.

%Select the correct days from the SWC and SOilT datasets to use in the antecedent calcs:
SWC_ForAnt=SWC_6h(time2_6h,z);   
SoilT_ForAnt=SoilT_6h(time2_6h,z);

%Calculation of the antecedent dataset for SWC used in Rr function
Nlag=4;
SWCant_Rr=zeros(Nt,Nz,Nlag);
for i=1:Nlag
    SWCant_temp=zeros(Nt,Nz);
    for j=1:28,  %28=7*4 because 7 days in the past and 4 six-hourly periods in each day.
        SWCant_temp=SWCant_temp+SWC_ForAnt((((Nt_beg-1)*4)+1)-((i-1)*28)-j:(Nt_end*4)-((i-1)*28)-j,1:Nz);
        %31 is the no. of days before day 1, where day 1 = 04/01/2007;
    end
    SWCant_Rr(:,:,i)=SWCant_temp/28;
end

%Calculation of the antecedent dataset for SWC used in Rm function:
Nlag=4;
SWCant_Rm=zeros(Nt,Nz,Nlag);
for i=1:Nlag
    SWCant_temp=zeros(Nt,Nz);
    for j=1:4,  %4 six-hourly periods in each day.
        SWCant_temp=SWCant_temp+SWC_ForAnt((((Nt_beg-1)*4)+1)-((i-1)*4)-j:((Nt_end)*4)-((i-1)*4)-j,1:Nz);
    end
    SWCant_Rm(:,:,i)=SWCant_temp/4;
end

%Calculation of the antecedent dataset for SoilT:
Nlag=4;
SoilTant=zeros(Nt,Nz,Nlag);
for i=1:Nlag
    SoilTant_temp=zeros(Nt,Nz);
    for j=1:4,  %4 six-hourly periods in each day.
        SoilTant_temp=SoilTant_temp+SoilT_ForAnt((((Nt_beg-1)*4)+1)-((i-1)*4)-j:(Nt_end*4)-((i-1)*4)-j,1:Nz);
    end
    SoilTant(:,:,i)=SoilTant_temp/4;
end

%Create vectors of the relative depth-dependent distributions of microbial biomass
%carbon (Mrel), root biomass (Rrel) and soil organic matter carbon (Crel). Then combine 
%into a matrix called Xrel.  The parameters are chosen based on data of microbial
%biomass carbon, root biomass and soil organic matter carbon.
if (Xrel_params==0)
    z = 1:1:100;      
    Mrel = gampdf(z,1.7,4.75);  %modelled distribution (by depth) of MBC (mg/cm3)
    Rrel=exp(z/-7);    %modelled distribution (by depth) of root biomass (mg/cm3)
    Crel=exp(-z/30);   %modelled distribution (by depth) of total C mass (mg/cm3
    Mrel(1,1:2)=repmat(Mrel(1,3),1,2);   %We replace the first 2 elements with the third.
    Rrel(1,1:2)=repmat(Rrel(1,3),1,2);
    Crel(1,1:2)=repmat(Crel(1,3),1,2);
    Mrel = Mrel/sum(Mrel);
    Rrel = Rrel/sum(Rrel);
    Crel = Crel/sum(Crel);
    Xrel=[Crel; Mrel; Rrel];
end

%Constants that are input to the soilCO2 function (ie don't change these):
dx = [6, 0.01];  %[dt dx] = time (hours) and depth (metres) steps.

% upper boundary conditions.
UPBC = atmCO2_IN;
%We don't have a lower boundary condition because we assume that dc/dz=0 at
%z=1 meter.

%Define depth-dependent distribution of d13Cr and d13m.  d13Cr and d13Cm are 
%the 13C signature of the CO2 respired by roots and microbes, respectively. 
%These are vectors that allow for a different signature at each depth/node, but 
%for Ryan et al. (in review) we assumed that that they are constant with depth. 
%We choose values for d13Cr and d13Cm below based on isotope data from PHACE site 
%which Elise Pendall provided. The data are for plot 11 (an elevated CO2 plot)
%which has the most similar SWC to plot 14 (an ambient CO2 and ambient warming plot 
%which was used in the computer experiment) compared to any of the other ECN plots.  
%We use an elevated CO2 plot because the difference between the d13Cr and d13Cm is 
%much larger than those for an ambient CO2 plot (which are around -24 for d13Cr and -19 for d13Cm).

z = 1:1:100; %here z refers to depth, so z=1 is 1cm depth, z=2 is 2cm depth, etc
d13Cr = -35*ones(1,101);  %Units are 0/00.  
d13Cm = -19*ones(1,101);  %no. of cols is 101 because we also include the soil's surface (0cm depth) 

%Store all the environmental and site data into an object called 'Data'.
Data.SWC = SWC;
Data.SoilT = SoilT;
Data.SWCantRm = SWCant_Rm;
Data.SWCantRr = SWCant_Rr;
Data.SoilTant = SoilTant;
Data.AtmPress = AtmPress;
Data.Gness = Gness;
Data.params_swp = params_swp;
Data.Xrel = Xrel;
Data.cinit = cinit;
Data.d13Cr = d13Cr;
Data.d13Cm = d13Cm;
Data.dx = dx;
Data.BC = UPBC;

%Write 'Data' to a mat file.  
if (write_file==1)
    save([num2str(site) 'data_6hourly.mat'],'Data')
end

if (Xrel_params==1)
    a0=12; a1=1.7; a2=4.75;
    Depth_Mobs=[2.5 10 22.5];
    Depth_all=[1:1:100];
    MBC_obs=[1.05 0.78 0];  %mg/cm3 of MBC for first three depths (DSS=1,2,3) for 2007.
    Mrel_obs=MBC_obs
    Mrel = a0*gampdf(Depth_Mobs,a1,a2)
    Mrel_all = a0*gampdf(Depth_all,a1,a2);
    Mrel_all(1,1:2)=repmat(Mrel_all(1,3),1,2);

    a0=18; a1=7;
    Depth_Robs=[2.5 10 22.5 37.5 60 87.5];
    Depth_all=[1:100];
    RootBiomass_obs=[12.9 3.29 0.868 0.253 0.0996 0.103];
    Rrel_obs=RootBiomass_obs
    Rrel=a0*exp(-Depth_Robs/a1)+0.0035
    Rrel_all=a0*exp(-Depth_all/a1);

    a0=25; a1=30;
    Depth_Cobs=[2.5 10 22.5];
    Depth_all=[1:100];
    soilC_obs=[22.26 19.5 11.29];  %mgC/cm3.
    Crel_obs=soilC_obs
    Crel=a0*exp(-Depth_Cobs/a1)+0.0035
    Crel_all=a0*exp(-Depth_all/a1);

    colr=['k' 'b' 'r' 'g' 'y'];

    f1=figure('Position',[10 40 1010 650]);           
    %This specifies the size of the figure when outputed on the screen;
    %figure('position',[x y w h]) x=width from LHS of screen to LHS of figure;
    %y=height from bottom of screen to bottom of figureq;
    %w = width of figure; h=height of figure;
    subplot(2,2,1)
    L1=plot(Depth_all,Mrel_all,'b-');
    hold on; 
    L2=plot(Depth_Mobs,Mrel_obs,'*b');
    hold on
    L3=plot(Depth_all,Rrel_all,'g-');
    hold on
    L4=plot(Depth_Robs,Rrel_obs,'*g');
    hold on
    L5=plot(Depth_all,Crel_all,'r-');
    hold on
    L6=plot(Depth_Cobs,Crel_obs,'*','color',[0.6,0.6,0.6]);
    hold on
    plot(Depth_Cobs,Crel_obs,'*r');
    hold on
    title('Vertical distn of MBC, root C and SOM C (mgC/cm^3)');
    legend([L5,L3,L1,L6],'SOM C','root C','MBC','data','Location',[0.3,0.8,0.1,0.1]);
    xlabel(['Depth (cm)'])
    ylabel(['C content (mg/cm^3)'])

    saveas(f1,['VerticalDistribution_ALL.jpg']) ;
    close all
end


    