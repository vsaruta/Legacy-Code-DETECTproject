clear all
close all

f1=figure('Position',[10 40 1010 650]);           
%This specifies the size of the figure when outputed on the screen;
%figure('position',[x y w h]) x=width from LHS of screen to LHS of figure;
%y=height from bottom of screen to bottom of figureq;
%w = width of figure; h=height of figure;
for s=1:3,
    ah(s)=subplot(3,1,s)
    if (s==1)
        ant=['']; titlename=['(a) \it{Ctrl} \rm{scenario, daily time-scale}'];
    elseif (s==2)
        ant=['_ant']; titlename=['(b) \it{Ctrl-ant} \rm{scenario, daily time-scale}'];
    elseif (s==3)
        titlename=['(c) \it{Ctrl-ant} \rm{scenario, sub-daily time-scale}']; 
    end
    days_IN=[91 135; 136 150; 151 160; 161 180; 181 216; 217 240; 241 273]-90;
        
    %Read in predicted CO2 from DETECT:
    if (s<=2)
        %Read in the temperature data for plot:
         load([pwd '/DETECT_Versions/SS_DETECT/SS_Inputs/PHACEdata_6hourly'  '.mat'])
%         load(['../Inputs/PHACEdata_6hourly.mat'])
        T_IN = Data.SoilT + 273.15;          %Nt * Nz
        T_3cm = T_IN(1:732,4);
        T_10cm = T_IN(1:732,11);
        T_20cm = T_IN(1:732,21);

        %CO2 preds:
%         load(['../Outputs/output_6hourly' num2str(ant) '.mat'])
        load([pwd '/DETECT_Versions/SS_DETECT/SS_Outputs/output_6hourly'  '.mat'])
        cctype=out.CO2type;
        cc=squeeze(sum(cctype(1:4,:,:),1));   %units: mgCO2 m-3  
        begT=1; endT=183*4;
        cc_3cm_temp=cc(4,2:733)'.*(8.314*T_3cm(begT:endT,1))./((44/1000)*101300);
        cc_10cm_temp=cc(11,2:733)'.*(8.314*T_10cm(begT:endT,1))./((44/1000)*101300);
        cc_20cm_temp=cc(21,2:733)'.*(8.314*T_20cm(begT:endT,1))./((44/1000)*101300);
        cc_3cm=mean(reshape(cc_3cm_temp',4,183)',2);
        cc_10cm=mean(reshape(cc_10cm_temp',4,183)',2);
        cc_20cm=mean(reshape(cc_20cm_temp',4,183)',2);        
        
        %Precip
        DayN1=[92:274]';  %DoY for 2008
        DayN2=[457:639]'; %2008 (DoYs 92 to 274) 
        precip_IN=csvread(['Data_for_plots/Precip_PHACE_fromKevin_2007to2013_11Mar2015.csv'],1,1);
        precip=precip_IN(DayN2,3);  %2008 precip
    
        %Read in the temperature data for plot 14ACN:
        T_IN = Data.SoilT + 273.15;          %Nt * Nz
        T_3cm = mean(reshape(T_IN(1:732,4),4,183))';
        T_10cm = mean(reshape(T_IN(1:732,11),4,183))';
        T_20cm = mean(reshape(T_IN(1:732,21),4,183))';
    
        %Read in soil CO2 data and for plot 14ACN and convert CO2 from ppm to mg CO2 m-3:
        cc_obs_ALL3plots=csvread(['Data_for_plots/Master Vaisala 2008 2009 flux calcs Ed ppm.csv'],1,0);
        cc_obs=cc_obs_ALL3plots(1:67,3:5);  %days 207 to 273;
        cc_obs_3cm = cc_obs(:,1); %days 207 to 273;
        cc_obs_10cm = cc_obs(:,2); %days 207 to 273;
        cc_obs_20cm = cc_obs(:,3); %days 207 to 273;
        DayN3=[207:273]'; 
    
    elseif(s==3)
        %Read in the temperature data for plot 14ACN:
        load([pwd '/DETECT_Versions/SS_DETECT/SS_Inputs/PHACEdata_6hourly'  '.mat'])
%         load(['../Inputs/PHACEdata_6hourly.mat'])
        T_IN = Data.SoilT + 273.15;          %Nt * Nz
        T_3cm = T_IN(1:732,4);
        T_10cm = T_IN(1:732,11);
        T_20cm = T_IN(1:732,21);
            
%         load(['../Outputs/output_6hourly' num2str(ant) '.mat'])
        load([pwd '/DETECT_Versions/SS_DETECT/SS_Outputs/output_6hourly'  '.mat'])
        cctype=out.CO2type;
        cc=squeeze(sum(cctype(1:4,:,:),1));   %units: mgCO2 m-3  
        begT=1; endT=183*4;
        cc_3cm=cc(4,2:733)'.*(8.314*T_3cm(begT:endT,1))./((44/1000)*101300);
        cc_10cm=cc(11,2:733)'.*(8.314*T_10cm(begT:endT,1))./((44/1000)*101300);
        cc_20cm=cc(21,2:733)'.*(8.314*T_20cm(begT:endT,1))./((44/1000)*101300);
    
        DayN1=[91.25:0.25:274]';  %DoY for 2008/2009 (DoYs 92 to 274)
        DayN2=[456.25:0.25:639]'; %2008  (days 457 to 639, or DoYs 92 to 274 for 2008)
        DayN2d=[457:639]';
        precip_IN=csvread(['Data_for_plots/Precip_PHACE_fromKevin_2007to2013_11Mar2015.csv'],1,1);
        precip=reshape(repmat(precip_IN(DayN2d,3)',4,1),1,4*length(DayN2d))';
    
        %Read in soil CO2 data and for plot 14ACN and convert CO2 from ppm to mg CO2 m-3:
        cc_obs_ALL3plots=csvread(['Data_for_plots/Master Vaisala 2008 2009 flux calcs Ed ppm hourly.csv'],1,0);
        cc_obs=cc_obs_ALL3plots(:,4:6);  %days 207 to 273;
        begT=206+(1/24); endT=273;
        cc_obs_3cm = cc_obs(:,1); %days 207 to 273;
        cc_obs_10cm = cc_obs(:,2); %days 207 to 273;
        cc_obs_20cm = cc_obs(:,3); %days 207 to 273;
        DayN3=[begT:(1/24):endT]'; 
    end
    if (s<=2)
        L7=bar(DayN1,precip'*100,0.5,'EdgeColor','b','LineWidth',1.5);
        hold on
        L1=plot(DayN1,cc_3cm,'-','color','k','LineWidth',0.75);
        hold on; 
        L2=plot(DayN1,cc_10cm,'-','color','g','LineWidth',0.75);
        hold on; 
        L3=plot(DayN1,cc_20cm,'-','color','r','LineWidth',0.75);
        hold on
        L4=plot(DayN3,cc_obs_3cm,'*','color','k','LineWidth',0.5);
        hold on; 
        L5=plot(DayN3,cc_obs_10cm,'*','color','g','LineWidth',0.5);
        hold on; 
        L6=plot(DayN3,cc_obs_20cm,'*','color','r','LineWidth',0.5);
        xlim([88 275])
    elseif (s==3)
        L7=bar(DayN1,precip'*100,0.5,'EdgeColor','b','LineWidth',1.5);
        hold on
        L1=plot(DayN1,cc_3cm,'-','color','k','LineWidth',0.5);
        hold on; 
        L2=plot(DayN1,cc_10cm,'-','color','g','LineWidth',0.5);
        hold on; 
        L3=plot(DayN1,cc_20cm,'-','color','r','LineWidth',0.5);
        hold on
        DayN3_3cm=DayN3(cc_obs_3cm>0); cc_obs_3cm=cc_obs_3cm(cc_obs_3cm>0);
        DayN3_10cm=DayN3(cc_obs_10cm>0); cc_obs_10cm=cc_obs_10cm(cc_obs_10cm>0);
        DayN3_20cm=DayN3(cc_obs_20cm>0); cc_obs_20cm=cc_obs_20cm(cc_obs_20cm>0);
        L4=plot(DayN3_3cm(1:2:length(cc_obs_3cm),1),cc_obs_3cm(1:2:length(cc_obs_3cm),1),'*','color','k','MarkerSize',2.25);
        hold on; 
        L5=plot(DayN3_10cm(1:2:length(cc_obs_10cm),1),cc_obs_10cm(1:2:length(cc_obs_10cm),1),'*','color','g','MarkerSize',2.25);
        hold on; 
        L6=plot(DayN3_20cm(1:2:length(cc_obs_20cm),1),cc_obs_20cm(1:2:length(cc_obs_20cm),1),'*','color','r','MarkerSize',2.25);
        xlim([200 274])
    end
    ylim([0 12500])
    if (s==1)
        set(gca, 'YAxisLocation','Left');
        set(gca, 'YTickLabel',{'0' '5000' '10000'});
        text(92,10500,[num2str(titlename)],'fontsize',12);
        text(276,5000,['50'],'fontsize',10);
        text(276,10000,['100'],'fontsize',10);
        legend([L1,L2,L3],'DETECT (3cm)','DETECT (10cm)','DETECT (20cm)','Location',[0.30,0.02,0.1,0.05]);
    elseif (s==2)
        ylabel(['Predicted CO_2 concentration (ppm)                      '],'fontsize',10) 
        text(92,10500,[num2str(titlename)],'fontsize',12);
        text(276,5000,['50'],'fontsize',10);
        text(276,10000,['100'],'fontsize',10);
        legend([L4,L5,L6],'Data (3cm)','Data (10cm)','Data (20cm)','Location',[0.50,0.02,0.1,0.05]);
    elseif (s==3)
        text(201,10500,[num2str(titlename)],'fontsize',11);
        text(274.4,5000,['50'],'fontsize',10);
        text(274.4,10000,['100'],'fontsize',10);
        legend([L7],'Precipitation','Location',[0.7,0.02,0.1,0.05]);
    end
    xlabel(['Day of year (Day 1 = 01/01/2008)'])
end

%find current position [x,y,width,height]
pos1 = get(ah(1),'Position');
pos2 = get(ah(2),'Position');
pos3 = get(ah(3),'Position');

%set width of 2nd, 3rd, etc... axes equal to first
pos2(3) = pos1(3);
pos3(3) = pos2(3);
set(ah(2),'Position',pos2);
set(ah(3),'Position',pos3);
set(ah(1),'XTickLabel','');

%set vertical position of 2nd, 3rd, etc... axes stacked on top of first.
pos2(2) = pos1(2) - pos2(4);
pos3(2) = pos2(2) - (1.5*pos3(4));
set(ah(2),'Position',pos2);
set(ah(3),'Position',pos3);

%insert tick marks etc... to right y-axis.
for i=1:3
    ax1 = axes('Position', get(ah(i),'Position'),'Color', 'none');
    set(ax1, 'YAxisLocation','Right','YTickLabel',{'0   ' '' '' '' '' '' ''});
    set(ax1, 'ticklength',[0 0])
    set(ax1, 'XTickLabel',{'' '' '' '' ''});
    if (i==2)
        set(ax1,'yaxislocation','right');
        ylabel(['Daily Precipitation (mm)                        '],'fontsize',10)
    end
end

saveas(f1,['SoilCO2_timeseries_PHACE.jpg']) ;
% close all


