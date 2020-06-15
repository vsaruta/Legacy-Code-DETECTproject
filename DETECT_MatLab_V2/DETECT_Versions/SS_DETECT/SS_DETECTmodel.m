function [out info] = DETECTmodel(Nt_run,SWC,SoilT,SWCant_Rm,SWCant_Rr,SoilTant,Gness,AtmPress,params_swp,Xrel,ic,d13Cr,d13Cm,dx,BC,params,params_lag,SS);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Based on soil CO2 production model developed by Fang & Moncrieff (1999)
% Agricultural and Forest Meteorology 95:225-236. But, simplify model to
% only deal with gas phase CO2 as is done in Hui & Luo (2004) Global
% Biogeochemical Cycles Vol 18:GB4029.
%
% INPUTS:
% Nt_sim = no. of time-steps to run the model for.  
% params, Xrel, SWC, SoilT, SWCant, SoilTant, Gness are used in SourceTerm.m.
% Xrel stores the distr. by depth of soil CO2, microb. C & root C (mgC/cm-3).
% SWC(t,z) = soil water content (m3/m3) at time t, depth z
% SoilT(t,z) = soil temp (C) at time t, depth z
% SWCant(t,z,lag) = soil water content (m3/m3) at time t - lag, depth z
% SoilTant(t,z) = soil temp (C) at time t - lag, depth z
% Gness(t,1) = Vegetation greeness at time t. 
% AtmPress(t,1) = atmostpheric pressure at the surface of PHACE.
% phig is the air filled porosity (m3/m3); phig100 is phig at -100cm H20.
% dx=[dt, dz]
% BC = [cppm, cdeep] (boundary conditions), where cppm = atmospheric [CO2]
% (or [CO2] at the soil surface) (ppm); cdeep = soil [CO2] at the bottom
% soil node (mg CO2 m-3).
% cinit(1,z) = initial conditions (soil [CO2], mg CO2 m-3) at time "zero"
% d13Cr = d13C of root-respired CO2, dim = 1xNz
% d13Cm = d13C of mirobial-respired CO2, dim = 1xNz
%
% OUTPUTS:
% out.CO2type =  is a matrix of [CO2] for every time (row), at every depth
% (col.) and CO2 type (12C-micr, 13C-micr, 12C-roots, 13C-roots).
% out.CO2 = the total [CO2], i.e. we sum across all four types.
% out.CO2fluxtype is a matrix of soil surface CO2 flux for every type (row)
% and time (col.).
% out.CO2flux = the total soil surface CO2 flux, i.e. we sum across all four types.
% out.d13C = total d13C of CO2.
% out.roots = soil surface CO2 flux from roots.

% Time step (dt) and depth increment (dz)
dt = dx(1);    % hours (hr)
dz = dx(2);    % meters (m)
% Number of time points (Nt), including initial time; and number of soil
% nodes (Nz), including surface layer (bottom boundary node is Nz+1)
[Nt Nz] = size(SWC);

%Create the phig dataset to be used when we compute the diffusion matrix 
%phig (air filled porosity) = phiT-SWC, where phiT (total soil porosity) = 1-(BD/PD)
%The value of SWC100 are estimated from fitting the water release curve function to site data of SWP.
BD=params_swp(1,:);
PD=2.52; 
SWC100=params_swp(3,:);
phiT=1-(BD./PD);  %air-filled porosity
phig=min(1,max(0,repmat(phiT,Nt,1) - SWC));  
phig100=repmat(phiT-SWC100,Nt,1);   

%Predict Diffusion at every depth & time using Moldrup's function. %Dgs0st = Diffusion 
%coefficient at standard atmosphere of temperature. T0=273.15K and pressure P0=101.3kPa.
Dg0st=repmat(0.0000139,Nt,Nz);  
T0 = 273.15;   
b=params_swp(2,:); 
P0=repmat(101.3,Nt,Nz);  %standard atmospheric pressure 
T = SoilT + 273.15;
P = repmat(AtmPress,1,Nz);  
logDg0= log(Dg0st)+(1.75*log(T/T0))+log(P0./P);
logDgs=logDg0+log((2*phig100.^3)+(0.04*phig100))+((2+(3./repmat(b,Nt,1))).*log(phig./phig100));
Dgs = exp(logDgs)*3600; % Convert Dgs above (m2 s-1) to Dgs below (m2 hr-1)
Dgs_13C = Dgs/1.0044;
Dgs_12C = Dgs;

%Upper boundary condition.  Note that we don't require a lower boundary conditions 
%because we assume dC/dz = 0 at the lower boundary (i.e. at 1m depth)  
catm = BC; % atm CO2 (ppm)
 
%We now need to compute the number of numerical time-steps needed to numerically solve the soil CO2 
%partial differential equation.  The formula given from Haberman, R. (1998). Elementary applied 
%partial differential equations (Third Edition), Prentice Hall Englewood Cliffs, NJ. (pg. 240) is
%Ndt = (dt*Dgs)/((dz^2)*0.5).  However, we need to do max(max(Dgs)) because Dgs is a matrix.  The 
%-0.05 %ensures that we have more numerical time-steps than the minimum, as a precaution.
Ndt = ceil(dt*max(max(Dgs))./((dz^2)*(0.5-0.05)));  

%Rall = Total microbial + root respiration at all depths and times. (mg C cm-3 hr-1):
%Note that Rall goes from depth z = 1cm to 100cm, so we add extra column of
%zeros (lines 91 to 93) to indicated that respiration = 0 at soil/atmosphere interface (z=0).
Rall = SS_SourceTerm(params,params_lag,Xrel,SWC,SoilT,SWCant_Rm,SWCant_Rr,SoilTant,Gness,SS);
R = [zeros(Nt,1) Rall.R];
Rr= [zeros(Nt,1) Rall.Rr];
Rm= [zeros(Nt,1) Rall.Rm];

% Convert to total source flux (mg CO2 m-3 hr-1) for each source and for each isotope 
% (13C and 12C): Current units for R are mg C cm-3 hr-1 per layer.  So, just need to (1) 
% convert from cm3 to m3 (mult. by 100^3); (2) convert from C to CO2 (mult. by 44/12).
S = R*(44/12)*(100^3);
Sr = Rr*(44/12)*(100^3);
Sm = Rm*(44/12)*(100^3);

%Partition the soil CO2 sources into the isotopic components:
%The lines below (115-120) are based on solving eqns (1) and (2) simultaneously for 
%[13C] (=S13root) and [12C] (=S12root):
%delta13Cr = [(Rr – Rstd)/Rstd] * 1000   ------------------- (1)
%where Rr = [13C] / [12C]
%[13C] + [12C] = Sr    ------------------------------------- (2)
%We solve for S13micr and S12micr in the same way.
%Note that the d13Cr and d13Cm refer to the isotopic signature of the
%micr. and root respiration at each depth, i.e. it's not hte isotopic
%signature of the surface respired CO2.  For the 2nd DETECT paper, we
%may wish to allow these d13Cr and d13Cm terms to vary with time (but 
%still keep fixed with depth).  
Rstd = 0.0112372;  % Rstd for 13C/12C (Pee Dee Belemite):
denomR = (1000+(1000+repmat(d13Cr,[Nt 1]))*Rstd);
denomM = (1000+(1000+repmat(d13Cm,[Nt 1]))*Rstd);
S13root = Sr.*(1000+repmat(d13Cr,[Nt 1])).*Rstd./denomR;
S12root = 1000*Sr./denomR;
S13micr = Sm.*(1000+repmat(d13Cm,[Nt 1])).*Rstd./denomM;
S12micr = 1000*Sm./denomM;
Stype(1,:,:)=S13root;
Stype(2,:,:)=S12root;
Stype(3,:,:)=S13micr;
Stype(4,:,:)=S12micr;
Stype(:,:,1)=0;

%Convert upper boundary condition (upbc) (ie atmospheric [CO2] or [CO2] at the soil
%surface) from ppm units to mg CO2 m-3 units.  Note that 101300 is the standard 
%pressure in Pa, R=8.134 is the gas constant, and the soil temp T is in Kelvin.
upbc = (catm*44/1000)*101300./(8.314*T(:,1)); 

%The Partial differential equation (equation 1 of Ryan et al., in review) is solved 
%four times for different types of soil CO2: (1) the d12C of soil CO2 from microbial 
%sources; (2) the d13C of soil CO2 from microbial sources; (3) the d12C of soil CO2 
%from root sources; (4) the d13C of soil CO2 from root sources.  This for
%loop computes the initial and boundary conditions for the four types of
%soil CO2.
cctype = ones(4,Nz+1,Nt+1)*-9999;
fluxtype=zeros(4,Nt);
for i=1:4,
    cctype(i,2:Nz,1) = (squeeze(Stype(i,1,2:Nz))'./S(1,2:Nz)).*ic(2:Nz);
    cctype(i,Nz+1,1) = cctype(i,Nz,1); %dC/dz = 0 at z=Nz means c(Nz+1,t)=c(Nz,t)
    cctype(i,1,1)=(squeeze(Stype(i,1,2))'./S(1,2)).*ic(1);
    cctype(i,1,2:(Nt+1)) = (squeeze(Stype(i,1,2))'./S(1,2)).*upbc;
end;

%Run main For Loop.  This computes the soil CO2 concentration for each of the four 
%types described above by numerically solving the partial differential equation.
ctype=zeros(4,Nz+1);
dt_old=dt;
dt = dt_old/Ndt;


%%% VVS COMMENT:
%%% Calculating Stype_Star, which is average over all depth and times
%%% for different sources (i)
for t=1:Nt_run,
    for z=1:Nz,
        for i=1:4,
            %Nz = 100 is 'N' of depth 'z'
            Stype_star(i,t) = sum(Stype(i,t,z))/(Nz);
            %Stype_star(i,t) = sum(Stype(i,t,z))/(Nz);
        end
    end
end

%%% VVS COMMENT:
% correction_argument DOES NOT HAVE ANY SCIENTIFIC SENSE JUST DEBUGGING COMPONENT.
% Set it to = 1 to use the formula from the paper
% catm OR upbc in other units the main driver of output for SS
correction_argument = 1;
% 0.3 = 293
% 0.35 = 586
% 0.4 = 879
% 0.334641638225256 ~ 496;

%%% VVS COMMENT:
% HERE SELECT WHICH 'FOOR LOOP' TO RUN FOR DEBUGGING PURPUSE
% 1 - At this point most simplified SS code (more then twice faster then
% others)
% 2 - Not simplified NSS efficient code with main calculation swapped
% to be SS (least intrusion code). Result is a little different from #1
% 3- Currently in development
which_to_run = 1; 



if which_to_run == 1,
%%% HERE STARTS SS LOOP
%%% VERY FAST LOOP ~2 seconds output is very similar to Loop 2
    
for t=1:Nt_run,
    % printing out current run N
    t
    for i =1:4,
        % If our Run is not first, we will recalculate previous time step
        % of cctype I wonder if this may cause issues, but not doing this
        % cause our results to be ~ 0.
        if t > 1,
            cctype(i,1,t)=(squeeze(cctype(i,2,t-1))'./sum(cctype(:,2,t-1))).*upbc(t-1);
            cctype(i,Nz+1,t)=cctype(i,Nz,t);
        end
%         ctype = squeeze(cctype(:,:,t));
        for z=2:Nz,
            % Dpending on time step and i set correct Dgs_T
             if ((i==1) || (i==3))
                Dgs_t= [Dgs_13C(t,1) Dgs_13C(t,:)];
                Dgs_tC = Dgs_13C;
             else
                Dgs_t= [Dgs_12C(t,1) Dgs_12C(t,:)];
                Dgs_tC = Dgs_12C;
             end
             % This SS solution code is based on Eq 11. of Ryan et al.
             % Stype_star for each source (i) and each time step (t) 
             % devided by Dgs_t(for each depth) multiplied by 'z-((z^2)/2)'
             % and add + upbc (which is 'catm' in different units)
              ctype(i,z) = (((Stype_star(i,t) / (Dgs_t(1,z))))*(z-((z^2)/2))) + 356; %(upbc(t,1)*correction_argument);
            % ctype(i,z) = (((Stype_star(i,t) / (Dgs_t(1,z))))*(0.1-((0.1^2)/2))) + (upbc(t,1)*correction_argument);
             %              ctype(i,z)
        end
    end
    cctype(:,:,t+1)=[ctype(1:4,1:Nz) ctype(1:4,Nz)];
    for i =1:4,
        MeanDgs = Dgs_tC(t,1);
        fluxtype(i,t) = (MeanDgs*(cctype(i,2,t)-cctype(i,1,t))/dz)/(44*3.6);
    end;
end;
%%% END OF MAIN LOOP


%%%% HERE STARTS ALTERNATIVE LOOP 2

elseif which_to_run == 2,
for t=1:Nt_run,
    t
    New_Dgs_13 = [Dgs_13C(t,1) Dgs_13C(t,:)]; 
    New_Dgs_12 = [Dgs_12C(t,1) Dgs_12C(t,:)];
    if t > 1,
        for i =1:4,
            cctype(i,1,t)=(squeeze(cctype(i,2,t-1))'./sum(cctype(:,2,t-1))).*upbc(t-1);
            cctype(i,Nz+1,t)=cctype(i,Nz,t);  %dC/dz = 0 at z=Nz means c(Nz+1,t)=c(Nz,t);
        end;
    end;
    ctype = squeeze(cctype(:,:,t));
    for tt=1:Ndt,
        ctype_old = ctype;
        for z=2:Nz,
            for i=1:4,
                if ((i==1) || (i==3))
                    Dgs_t= New_Dgs_13; 
                else
                    Dgs_t= New_Dgs_12;
                end
                % Numerical approximation of CO2 concentrations based on the PDE eqn in
                % Fang & Moncrieff. Numberical sol'n based on the nonhomogeneous problem 
                % in Haberman (1998) Elementary Applied Partial Differential Equations, 
                % Page 240 (section 6.3.7), aka Forward Euler method:
                if (z==Nz)
%                     ctype(i,z) = ctype_old(i,z) +...
%                         dt*((Dgs_t(1,z)/(dz*dz))*(ctype_old(i,z-1)-ctype_old(i,z)) +...
%                         Stype(i,t,z));
                      ctype(i,z) = (((Stype_star(i,t) / (Dgs_t(1,z))))*(z-((z^2)/2))) + (upbc(t,1)*correction_argument);
                else
%                     ctype(i,z) = ctype_old(i,z) +...
%                         dt*((Dgs_t(1,z)/(dz*dz))*(ctype_old(i,z+1)-2*ctype_old(i,z)+ctype_old(i,z-1)) +...
%                         ((Dgs_t(1,z+1)-Dgs_t(1,z-1))*(ctype_old(i,z+1)-ctype_old(i,z-1)))/(4*dz^2) +...
%                         Stype(i,t,z));     
                      ctype(i,z) = (((Stype_star(i,t) / (Dgs_t(1,z))))*(z-((z^2)/2))) + (upbc(t,1)*correction_argument);
                end
            end
        end
    end
    cctype(:,:,t+1)=[ctype(1:4,1:Nz) ctype(1:4,Nz)];
    for i=1:4,
        if ((i==1) || (i==3))
            MeanDgs=Dgs_13C(t,1);
        else
            MeanDgs=Dgs_12C(t,1);
        end
        fluxtype(i,t) = (MeanDgs*(cctype(i,2,t)-cctype(i,1,t))/dz)/(44*3.6);  
        %dividing by (44*3.6) converts units of soil CO2***flux from mg CO2 m-2 hr-1 to
        %umol CO2 m-2 s-1.
    end;
end;
%%%% END OF LOOP 2

%%%% HERE STARTS LOOP 3 IN DEVELOPMENT
%%%% Higly simplified buT CURRENTLY NO OUTPUT (NOT WORKING)
elseif which_to_run == 3,
disp('test')
i = 1
for t = 1:732,
    t
    for z = 1:100,
        % Dgs_t= [Dgs_13C(t,1) Dgs_13C(t,:)];
        Dgs_t = Dgs_13C(t,z);
        MeanDgs=Dgs_13C(t,1);
        %check_z_arg = z-((z^2)/2)
        ctype(i,z,t) = (((Stype_star(i,t) / (Dgs_t)))*(z-((z^2)/2))) + upbc(t,1);
        fluxtype(i,t) = (MeanDgs*(ctype(i,2,t)-ctype(i,1,t))/dz)/(44*3.6);
    end
end
%%% END OF MAIN LOOP 3

%%%% HERE STARTS LOOP 4 IN DEVELOPMENT
%%%% Expremental loop for tests
elseif which_to_run == 4,
    
%%%% Below is old version of SS

for t=1:Nt_run,
    for z=2:Nz,
        for i=1:4,
Stype_star(t) = sum(Stype(i,t,z+1))/(100);
        end
    end
end
for t=1:Nt_run,
    t
    New_Dgs_13C =[Dgs_13C(t,1) Dgs_13C(t,:)];
    New_Dgs_12C =[Dgs_12C(t,1) Dgs_12C(t,:)];
    if t > 1,
        for i =1:4,
            cctype(i,1,t)=(squeeze(cctype(i,2,t-1))'./sum(cctype(:,2,t-1))).*upbc(t-1);
            cctype(i,Nz+1,t)=cctype(i,Nz,t);  %dC/dz = 0 at z=Nz means c(Nz+1,t)=c(Nz,t);
            
        end;
    end;
    ctype = squeeze(cctype(:,:,t));
    for tt=1:Ndt,
        ctype_old = ctype;
        for z=2:Nz,
            for i=1:4,
                if ((i==1) || (i==3))
                    Dgs_t=New_Dgs_13C;
                else
                    Dgs_t=New_Dgs_12C;
                end
                % Numerical approximation of CO2 concentrations based on the PDE eqn in
                % Fang & Moncrieff. Numberical sol'n based on the nonhomogeneous problem 
                % in Haberman (1998) Elementary Applied Partial Differential Equations, 
                % Page 240 (section 6.3.7), aka Forward Euler method:
                
%                 Stype_star(t) = sum(Stype(i,t,z+1))/(100);
                if (z==Nz)
                    ctype(i,z) = (Stype_star(t)/Dgs_t(1,z))*(z-((z^2)/2))+catm(t,1);
%ctype(i,z) = ctype_old(i,z) +  dt*((Dgs_t(1,z)/(dz*dz))*(ctype_old(i,z-1)-ctype_old(i,z)) + Stype(i,t,z))
                else
                    ctype(i,z) = (Stype_star(t)/Dgs_t(1,z))*(z-((z^2)/2))+catm(t,1);                    
                end
            end
        end
    end
    cctype(:,:,t+1)=[ctype(1:4,1:Nz) ctype(1:4,Nz)];
    for i=1:4,
        if ((i==1) || (i==3))
            MeanDgs=Dgs_13C(t,1);
        else
            MeanDgs=Dgs_12C(t,1);
        end
        fluxtype(i,t) = (MeanDgs*(cctype(i,2,t)-cctype(i,1,t))/dz)/(44*3.6);  
        %dividing by (44*3.6) converts units of soil CO2 from mg CO2 m-2 hr-1 to
        %umol CO2 m-2 s-1.
    end;
end;  
    
%%% END OF MAIN LOOP 4


%%%% HERE STARTS LOOP 5 IN DEVELOPMENT
%%%% THIS IS NSS CODE + TRANSFORMATION OF NSS TO FAKE-SS (from Kim).
%%%% IT IS EDITED AND REARRANGED
elseif which_to_run == 5,
    for t=1:Nt_run,
    t
    New_Dgs_13 = [Dgs_13C(t,1) Dgs_13C(t,:)]; 
    New_Dgs_12 = [Dgs_12C(t,1) Dgs_12C(t,:)];
    if t > 1,
        for i =1:4,
            cctype(i,1,t)=(squeeze(cctype(i,2,t-1))'./sum(cctype(:,2,t-1))).*upbc(t-1);
            cctype(i,Nz+1,t)=cctype(i,Nz,t);  %dC/dz = 0 at z=Nz means c(Nz+1,t)=c(Nz,t);
        end;
    end;
    ctype = squeeze(cctype(:,:,t));
    for tt=1:Ndt,
        ctype_old = ctype;
        for z=2:Nz,
            for i=1:4,
                if ((i==1) || (i==3))
                    Dgs_t= New_Dgs_13; 
                else
                    Dgs_t= New_Dgs_12;
                end
                % Numerical approximation of CO2 concentrations based on the PDE eqn in
                % Fang & Moncrieff. Numberical sol'n based on the nonhomogeneous problem 
                % in Haberman (1998) Elementary Applied Partial Differential Equations, 
                % Page 240 (section 6.3.7), aka Forward Euler method:
                if (z==Nz)
                    ctype(i,z) = ctype_old(i,z) +...
                        dt*((Dgs_t(1,z)/(dz*dz))*(ctype_old(i,z-1)-ctype_old(i,z)) +...
                        Stype(i,t,z));
                else
                    ctype(i,z) = ctype_old(i,z) +...
                        dt*((Dgs_t(1,z)/(dz*dz))*(ctype_old(i,z+1)-2*ctype_old(i,z)+ctype_old(i,z-1)) +...
                        ((Dgs_t(1,z+1)-Dgs_t(1,z-1))*(ctype_old(i,z+1)-ctype_old(i,z-1)))/(4*dz^2) +...
                        Stype(i,t,z));                    
                end
            end
        end
    end
    cctype(:,:,t+1)=[ctype(1:4,1:Nz) ctype(1:4,Nz)];
    for i=1:4,
        if ((i==1) || (i==3))
            MeanDgs=Dgs_13C(t,1);
        else
            MeanDgs=Dgs_12C(t,1);
        end
        fluxtype(i,t) = (MeanDgs*(cctype(i,2,t)-cctype(i,1,t))/dz)/(44*3.6);  
        %dividing by (44*3.6) converts units of soil CO2 from mg CO2 m-2 hr-1 to
        %umol CO2 m-2 s-1.
    end;
end;

        phi_SS=squeeze(sum(Stype([1:4],:,:),1));
        Smicrobes_IN=squeeze(sum(Stype([3 4],:,:),1)); 
        Sm = mean(Smicrobes_IN(1:732,2:100),2);  %Source term averaged over all depths.
        Sroots_IN=squeeze(sum(Stype([1 2],:,:),1)); 
        Sr = mean(Sroots_IN(1:732,2:100),2);  %Source term averaged over all depths.
        Rstd=0.0112372; d13Catm=-19;  %from Elise.
        
        d13Cr=-35; d13Cm=-19; 
        Ar=1000/(((d13Cr+1000)*Rstd)+1000);
        Br=((d13Cr+1000)*Rstd)/(((d13Cr+1000)*Rstd)+1000);
        Am=1000/(((d13Cm+1000)*Rstd)+1000);
        Bm=((d13Cm+1000)*Rstd)/(((d13Cm+1000)*Rstd)+1000);
        
        S12 = (Ar*Sr) + (Am*Sm);
        S13 = (Br*Sr) + (Bm*Sm);
        
        Ds12_1cm = Dgs_12C(1:732,1);
        Ds13_1cm = Dgs_13C(1:732,1);
        
        z=0.01;
        dz=0.01;
        
        %Compute Catm12 and Catm13
        %out.CO2type == cctype
        Catm12_IN_v1=squeeze(sum(cctype([2 4],1,2:733),1)); 
        Catm13_IN_v1=squeeze(sum(cctype([1 3],1,2:733),1)); 
        %
        Catm_IN=squeeze(sum(cctype([1:4],1,2:733),1));   %mg CO2 m-3;     
        Catm12_IN=(1000*Catm_IN)/(((d13Catm+1000)*Rstd)+1000);
        Catm13_IN=((Catm_IN*(d13Catm+1000)*Rstd))/(((d13Catm+1000)*Rstd)+1000);
        %We now compute teh steady-state solution for soil CO2 at a depth
        %of 1cm:
        Catm13=Catm13_IN(1:732,1);
        Catm12=Catm12_IN(1:732,1);
        %
        C12_1cm = ((S12./Ds12_1cm).*z.*(1-(z/2))) + Catm12;
        C13_1cm = ((S13./Ds13_1cm).*z.*(1-(z/2))) + Catm13;
        Flux_C13=Ds13_1cm.*((C13_1cm-Catm13)/dz); 
        Flux_C12=Ds12_1cm.*((C12_1cm-Catm12)/dz);           
        FluxSS_temp=(Flux_C13+Flux_C12)/(44*3.6);  %units of FluxSS_temp: umol CO2 m-2 s-1
        FluxSS=mean(reshape(FluxSS_temp(1:732,1),4,183),1)'
        
        
%         cctype(:,:,t+1)=[ctype(1:4,1:Nz) ctype(1:4,Nz)];
%         for i=1:4,
%             if ((i==1) || (i==3))
%                 MeanDgs=Dgs_13C(t,1);
%             else
%                 MeanDgs=Dgs_12C(t,1);
%             end
%         fluxtype(i,t) = (MeanDgs*(cctype(i,2,t)-cctype(i,1,t))/dz)/(44*3.6);  
%         fluxtype = FluxSS;

%%% END OF MAIN LOOP 5
end; %AND of which to RUN!



%Store all the important outputs to 'out'.
out.CO2type = cctype;
out.CO2 = squeeze(sum(cctype,1))';
out.CO2fluxtype = fluxtype;
out.CO2flux = sum(fluxtype,1)';
out.d13C = real(squeeze((((sum(fluxtype([1 3],:),1)./sum(fluxtype([2 4],:),1))./Rstd)-1)*1000)');
out.roots = real(squeeze(sum(fluxtype(1:2,:),1)./sum(fluxtype,1))');

%Store all the less important model output as 'info'. 
info.Dgs=Dgs; info.phig=phig; info.R=R; info.S=S; 
info.Rm=Rm; info.Rr=Rr; info.Stype=Stype;