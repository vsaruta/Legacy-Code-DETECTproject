clear all
close all

Ns=1; 
%Set to 1 if just plotting Rsoil predicted from the DETECT model WITHOUT the antecedent parameters.
%Set to 2 if plotting Rsoil predicted from the DETECT model WITH and WITHOUT the antecedent parameters.

%units from output files:
%Dgs: m2/hr
%S: mgCO2/m3/hr
%cc, cctype: mgCO2/m3
%CO2flux, CO2fluxtype: umolCO2/m3/s
%umolCO2/m2/s = (12*24*3600)/10^6 * gC/m2/day.
%C12=(S12/Ds12)*z*(1-z/2) (units equal on both sides).

%Read in daily ecosystem respiration (Reco) data.  We select the midday (11am, 12pm or 1pm) 
%Reco ecosystem chamber measurements and then convert them to daily Reco using the rule the 
%Elise Pendall used in her 2014 paper (see supplemental) which had an R2 of 0.93.  
DayN_RecoObs=[789 860 871 885 888 899 920 936 955 979] - (365*2);
RecoObs=[0.37 0.96 1.77 5.65 5.16 6.35 2.49 1.85 7.03 3.89]; %units: gC/m2/day
DayN_RmObs=[170 184 198 213 226 235 253 269];
RmObs=[1.48 1.22 0.52 0.46 2.54 1.61 1.31 1.01]; %units: gC/m2/day

%Define time-steps.
DayN1=[92:274]';  %DoY for 2008/2009
DayN2=[457:639]'; %2008

%Read in precipitation data:
precip_IN=csvread(['Data_for_plots/Precip_PHACE_fromKevin_2007to2013_11Mar2015.csv'],1,1);
precip=repmat(precip_IN(DayN2,3),1,4);

f1=figure('Position',[10 40 1010 650]);           
%This specifies the size of the figure when outputed on the screen;
%figure('position',[x y w h]) x=width from LHS of screen to LHS of figure;
%y=height from bottom of screen to bottom of figureq;
%w = width of figure; h=height of figure;
for s=1:Ns
    ah(s)=subplot(2,1,s)
    if (s==1)
        ant=['']; titlename=['(a) DETECT model without antecedent parameterization (\it{Ctrl} \rm{scenario)}'];
    elseif (s==2)
        ant=['_ant']; titlename=['(b) DETECT model with antecedent parameterization (\it{Ctrl-ant} \rm{scenario)}'];
    end
    
   %Load data:
%    load(['../Outputs/output_6hourly' num2str(ant) '.mat'])
   load([pwd '/DETECT_Versions/SS_DETECT/SS_Outputs/output_6hourly'  '.mat'])
   R_IN=out.CO2flux';
   Rr_IN=sum(out.CO2fluxtype(1:2,:));
   Rm_IN=sum(out.CO2fluxtype(3:4,:));
   R_DETECT=mean(reshape(R_IN(1,1:732),4,183),1);
   Rr_DETECT=mean(reshape(Rr_IN(1,1:732),4,183),1);
   Rm_DETECT=mean(reshape(Rm_IN(1,1:732),4,183),1);
        
   L1=bar(DayN1(1:183,1),precip(1:183,1)/5,0.5,'EdgeColor','b');
   hold on
   L2=plot(DayN1(1:183,1)',Rr_DETECT(1,1:183)*1.0368,'-','color','g','LineWidth',0.5); %*1.0368 converts to gC/m2/day.
   hold on; 
   L3=plot(DayN1(1:183,1)',Rm_DETECT(1,1:183)*1.0368,'-','color','b','LineWidth',0.5); %*1.0368 converts to gC/m2/day.
   hold on
   L4=plot(DayN1(1:183,1)',R_DETECT(:,1:183)*1.0368,'-','color','r','LineWidth',1); %*1.0368 converts to gC/m2/day.
   hold on
   L5=plot(DayN_RecoObs,RecoObs,'o','color','r','MarkerSize',8,'LineWidth',1) %units already in gC/m2/day.
   hold on
   L6=plot(DayN_RmObs,RmObs,'s','color','b','MarkerSize',8,'LineWidth',1) %units already in gC/m2/day.
   hold on
   L7=bar(279.5,sum(R_DETECT(:,1:183)*1.0368)/60,5,'FaceColor',[0.4 0.4 0.4],'EdgeColor',[0.4 0.4 0.4]);
   xlim([88 290])
   ylim([0 20])
   set(gca,'YTick',[0:5:15],'fontsize',10);
   set(gca,'XTick',[100:20:260],'fontsize',10);
   if (Ns==1)
       legend([L4,L3,L2,L5,L6,L1],'Total R_{soil} from DETECT','Root contribution','Microbial contribution','Observed daily R_{eco}','Observed daily R_m','Precipitation','Location',[0.475,0.3,0.1,0.15]);
       ylabel(['Daily respiration of CO_2 (gC m^{-2} day^{-1})'])
       text(275,-1.5,['Annual R_{soil} (gCm^{-2})'],'fontsize',9);
   elseif (Ns==2)
       if (s==1)
           legend([L4,L3,L2],'Total R_{soil} from DETECT','Root contribution','Microbial contribution','Location',[0.35,0.05,0.1,0.05]);
           ylabel(['Daily respiration of CO_2 (gC m^{-2} day^{-1})                                                     '])
       elseif (s==2)
           text(275,-1.5,['Annual R_{soil} (gCm^{-2})'],'fontsize',9);
           legend([L5,L6,L1],'Observed daily R_{eco}','Observed daily R_m','Precipitation','Location',[0.58,0.05,0.1,0.05]);  
       end
   end
   xlabel(['Day of year (Day 1 = 01/01/2008)'])
   text(92,16.5,[num2str(titlename)],'fontsize',12);
   text(277,(sum(R_DETECT*1.0368)/60)+1,[num2str(round(sum(R_DETECT*1.0368)))],'fontsize',8.5);
   text(291.5,0,['0'])
   text(291.5,5,['25'])
   text(291.5,10,['50'])
   text(291.5,15,['75'])
end
if (Ns==2)
    %find current position [x,y,width,height]
    pos1 = get(ah(1),'Position');
    pos2 = get(ah(2),'Position');
    
    %set width of second axes equal to first
    pos2(3) = pos1(3);
    set(ah(2),'Position',pos2);
    set(ah(1),'XTickLabel','');
    pos2(2) = pos1(2) - pos2(4);
    set(ah(2),'Position',pos2);
    
    %set width of second axes equal to first
    pos2(3) = pos1(3);
    set(ah(2),'Position',pos2);
    set(ah(1),'XTickLabel','');
    pos2(2) = pos1(2) - pos2(4);
    set(ah(2),'Position',pos2);
end

for i=1:Ns
    ax1 = axes('Position', get(ah(i),'Position'),'Color', 'none');
    set(ax1, 'YAxisLocation','Right','YTickLabel',{'    ' '    ' '    ' '    ' '    '});
    set(ax1, 'ticklength',[0 0])
    set(ax1, 'XTickLabel',{'' '' '' '' ''});
    if (i==2)
        set(ax1,'yaxislocation','right');
        ylabel(['                                                   Daily Precipitation (mm)'])
    end
end

if (Ns==1)
    saveas(f1,['Plots/SS_Rsoil_timeseries_PHACEsite.jpg'])
else
    saveas(f1,['Rsoil_timeseries_withAnt_PHACEsite.jpg'])
end
% close all