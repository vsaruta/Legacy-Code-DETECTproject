
function simpleGUI
clear all
close all
    padding_x = 15;
    padding_y = 15;
    first_panel_shift_y = 50;
    colmn_1_size = 300;
    colmn_2_size = 250;
    colmn_2_start = 220;
    colmn_3_start = 550;
    colmn_step = 55;
    pos_step = 25;
    x_var_col_size = 50;
    y_var_row_size = 20;
    hFig = figure('Visible','off', 'Menu','none', 'Name','DETECT models', 'Resize','off', 'Position',[100 100 (padding_x*2)+colmn_3_start+(colmn_step*5)+50 225+first_panel_shift_y]);
    movegui(hFig,'center')          %# Move the GUI to the center of the screen
%                        fig = uifigure;
%                        ax = uiaxes(hFig);
%                        x = linspace(-pi,pi,50);
% f = figure;
% ax = axes(hFig);
% ax.Units = 'pixels';
% ax.Position = [75 75 325 280]
% c = uicontrol;
% c.String = 'Plot Data';
% c.Callback = @plotButtonPushed;
    
    hBtnGrp = uibuttongroup('Position',[0 0 0.25 1], 'Units','Normalized');
    title_01 = uicontrol('Style','text', 'Position',[padding_x first_panel_shift_y+185 115 25], 'String','Choose Model to run:');
    
    
    
    
%     uicontrol('Style','frame', 'Parent',hBtnGrp, 'HandleVisibility','off', 'Position',[5 5 5 5], 'String','MatLab DETECT Original NSS', 'Tag','Model_1')

    uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off', 'Position',[padding_x first_panel_shift_y+150 colmn_1_size 30], 'String','MatLab DETECT Original NSS', 'Tag','Model_1')
    uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off', 'Position',[padding_x first_panel_shift_y+120 colmn_1_size 30], 'String','MatLab DETECT Efficient NSS', 'Tag','Model_2')
    uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off', 'Position',[padding_x  first_panel_shift_y+90 colmn_1_size 30], 'String','MatLab DETECT FAKE SS', 'Tag','Model_3')
    uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off', 'Position',[padding_x  first_panel_shift_y+60 colmn_1_size 30], 'String','MatLab DETECT SS (WiP)', 'Tag','Model_4')
    uicontrol('Style','Radio', 'Parent',hBtnGrp, 'HandleVisibility','off', 'Position',[padding_x  first_panel_shift_y+30 colmn_1_size 30], 'String','Rscript DETECT Efficient NSS', 'Tag','Model_5')
hBtnGrp2 = uibuttongroup('Position',[0.25 0 0.3 1], 'Units','Normalized');
title_02 = uicontrol('Style','text', 'Position',[padding_x+colmn_2_start first_panel_shift_y+185 115 25], 'String','Models Description:');    
    disc_M1 =  uicontrol('Style','text', 'Position',[padding_x+colmn_2_start first_panel_shift_y+150 colmn_2_size 25], 'String','- basic original model withou changes');
    disc_M2 =  uicontrol('Style','text', 'Position',[padding_x+colmn_2_start first_panel_shift_y+120 colmn_2_size 25], 'String','- 60 times more efficient model version');
    disc_M3 =  uicontrol('Style','text', 'Position',[padding_x+colmn_2_start first_panel_shift_y+90 colmn_2_size 25], 'String','- keep SWC, Gness, and Tsoil constant');
    disc_M4 =  uicontrol('Style','text', 'Position',[padding_x+colmn_2_start first_panel_shift_y+60 colmn_2_size 25], 'String','- simplified version with stable state solution');
    disc_M4 =  uicontrol('Style','text', 'Position',[padding_x+colmn_2_start first_panel_shift_y+30 colmn_2_size 25], 'String','- same as DETECT Efficient but in R code');
   hEdit3 = uicontrol('Style','text', 'Position',[padding_x+colmn_2_start first_panel_shift_y colmn_2_size 20], 'String','Status: Ready'); 
    uicontrol('Style','pushbutton', 'String','Run Model', 'Position',[padding_x first_panel_shift_y 200 25], 'Callback',{@button_callback})
%     hEdit3 = uicontrol('Style','text', 'Position',[padding_x first_panel_shift_y-35 200 20], 'String','Status: Ready');
    
%     [colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size]
hBtnGrp2 = uibuttongroup('Position',[0.55 0 0.45 1], 'Units','Normalized');
title_03 = uicontrol('Style','text', 'Position',[colmn_3_start first_panel_shift_y+185 115 25], 'String','Parameters Values:');
y_row_position = 0;
x_col_position = 0;
    param_1_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*0) x_var_col_size y_var_row_size], 'String','Rstar');
x_col_position = 1;
    param_1 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*0) x_var_col_size y_var_row_size], 'String','111.54');
x_col_position = 2;
    param_2_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*0) x_var_col_size y_var_row_size], 'String','RrBase');
x_col_position = 3;
    param_2 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*0) x_var_col_size y_var_row_size], 'String','0.00006');
x_col_position = 4;
    param_3_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*0) x_var_col_size y_var_row_size], 'String','Rr_alpha1');
x_col_position = 5;
    param_3 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*0) x_var_col_size y_var_row_size], 'String','11.65');

    
y_row_position = 1;
x_col_position = 0;
    param_4_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','Rr_alpha2');
x_col_position = 1;
    param_4 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','20.72');
x_col_position = 2;
    param_5_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','Rr_alpha3');
x_col_position = 3;
    param_5 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','-164.16');
x_col_position = 4;
    param_6_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','Cstar');
x_col_position = 5;
    param_6 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','711.6094');

    
y_row_position = 2;
x_col_position = 0;
    param_7_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','Mstar');
x_col_position = 1;
    param_7 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','12.2568');
x_col_position = 2;     
    param_8_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','RmBase');
x_col_position = 3;
    param_8 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','0.0015');
x_col_position = 4;    
    param_9_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','alpha1');
x_col_position = 5;
    param_9 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','14.05');

    
y_row_position = 3;
x_col_position = 0;
    param_10_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','alpha2');
x_col_position = 1;
    param_10 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','11.05');
x_col_position = 2;
    param_11_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','alpha3');
x_col_position = 3;
    param_11 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','-87.55');
x_col_position = 4;
    param_12_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','Km');
x_col_position = 5;
    param_12 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','0.00001');


y_row_position = 4;
x_col_position = 0;
    param_13_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','CUE');
x_col_position = 1;
    param_13 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','0.8');
x_col_position = 2;  
    param_14_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','solfrac');
x_col_position = 3;
    param_14 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','0.0004');
x_col_position = 4;
    param_15_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','Dliq');
x_col_position = 5;
    param_15 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','3.17');


y_row_position = 5;
x_col_position = 0;
    param_16_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','Eo_Rm');
x_col_position = 1;
    param_16 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','324.6');
x_col_position = 2;
    param_17_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','Eo_Rr');
x_col_position = 3;
    param_17 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','324.6');
x_col_position = 4;
    param_18_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','To');
x_col_position = 5;
    param_18 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','227.5');


y_row_position = 6;
x_col_position = 0;
    param_19_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','alpha4');
x_col_position = 1;
    param_19 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','-4.67');
x_col_position = 2;
    param_20_lable = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','tau');
x_col_position = 3;
    param_20 = uicontrol('Style','edit', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size y_var_row_size], 'String','10');
x_col_position = 4;  
    play_sound_leb = uicontrol('Style','text', 'Position',[colmn_3_start+(colmn_step*x_col_position) first_panel_shift_y+padding_y+(pos_step*y_row_position) x_var_col_size+10 y_var_row_size], 'String','Play Sound');
x_col_position = 5;    
    play_sound = uicontrol('Style','checkbox', 'Position',[colmn_3_start+(colmn_step*x_col_position)+15 first_panel_shift_y+padding_y+(pos_step*y_row_position)+3 x_var_col_size y_var_row_size]);

    hEdit4 = uicontrol('Style','text', 'Position',[colmn_3_start+235 first_panel_shift_y+185 115 25], 'String','Time to run: 0');
    hEdit5 = uicontrol('Style','text', 'Position',[colmn_3_start+150 first_panel_shift_y+185 115 25], 'String','Time to run:');
%     [padding_x+200 first_panel_shift_y 100 20]

% HERE STARTS SECOND PANNEL

uicontrol('Style','pushbutton', 'String','Plot Results', 'Position',[padding_x 15 200 25], 'Callback',{@button_callback2})

    set(hFig, 'Visible','on')        %# Make the GUI visible
    set(hEdit4, 'Visible','off')
    set(hEdit5, 'Visible','off')
    running_antecedent_vers=0;  % Including (=1) or not (=0) the antecedent variables in the DETECT model.
% runmodel=1;                 % Whether to run the model (=1) or not(=0).
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

% save([pwd '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/output_6hourly' num2str(ant) '.mat'],'out')
Data = [];
load([pwd '/DETECT_Versions/Efficient_DETECT/Efficient_Inputs/PHACEdata_6hourly.mat'])  %Control treatment plot from PHACE site for 2008
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


% params(1) = 111.54;  % Rstar = total root biomass (mgC) mass in a one cm2 column of soil 1m deep. 
params(1) = str2double(get(param_1, 'String'));

% params(2) = 0.00006; % RrBase, mass-specific root respiration base rate at a tau degC (mgC cm-3 hr-1).
params(2) = str2double(get(param_2, 'String'));

% params(3) = 11.65;    % alpha1 for Rr; SWC param.   
params(3) = str2double(get(param_3, 'String'));

% params(4) = 20.72;   % alpha2 for Rr; The main effect of SWCant on the log-scaled Rrbase
params(4) = str2double(get(param_4, 'String'));

% params(5) = -164.16; % alpha3 for Rr; The main effect of SWC*SWCant on the log-scaled Rrbase
params(5) = str2double(get(param_5, 'String'));


%Parameters used by the Rm submodel.
% params(6) = 711.6094; % Cstar = total soil C (mgC) mass in a one cm2 column of soil 1m deep. 
params(6) = str2double(get(param_6, 'String'));

% params(7) = 12.2568;  % Mstar = total MBC (mgC) mass in a one cm2 column of soil 1m deep. 
params(7) = str2double(get(param_7, 'String'));

% params(8) = 0.0015;   % RmBase; (mg C cm-3 hr-1); based on value used in Davidson et al. (2012)
params(8) = str2double(get(param_8, 'String'));

% params(9) = 14.05;    % alpha1 (mg C cm-3 hr-1); based on value used in Ryan et al. (2015) 
params(9) = str2double(get(param_9, 'String'));

% params(10) = 11.05;    % alpha2; The main effect of SWCant on the log-scaled Rmbase;  
params(10) = str2double(get(param_10, 'String'));

% params(11) = -87.55;   % alpha3; The main effect of SWC*SWCant on the log-scaled Rmbase; 
params(11) = str2double(get(param_11, 'String'));

% params(12) = 0.00001; % Km
params(12) = str2double(get(param_12, 'String'));

% params(13) = 0.8;     % CUE; Carbon Use Efficiency (mgC/mgC).  Taken from Tucker et al. (2012) 
params(13) = str2double(get(param_13, 'String'));

% params(14) = 0.0004;  % solfrac; Fraction of substrate that Sx that is soluble (called 'p' in DAMM paper). 
params(14) = str2double(get(param_14, 'String'));

% params(15) = 3.17;    % Dliq; Diffusivity of substrate Sx that is liquid.  Dimensionless.
params(15) = str2double(get(param_15, 'String'));

% params(16) = 324.6;    % Eo_Rm; Eo term used in Rm submodel (Kelvin). Use value from Ryan et al. (2015)
params(16) = str2double(get(param_16, 'String'));

% params(17) = 324.6;   % Eo_Rr; Eo term used in Rr submodel (Kelvin); Use value from Ryan et al. (2015)
params(17) = str2double(get(param_17, 'String'));

% params(18) = 227.5;   % To; Temp. parameter used in the Lloyd & Taylor temp. (Kelvin)
params(18) = str2double(get(param_18, 'String'));

% params(19) = -4.67;    % alpha4; The main effect of SoilTant on the log-scaled Rmbase;
params(19) = str2double(get(param_19, 'String'));

% params(20) = 10;      % tau; Base temperature in Kelvin.  We use 10 degC.
params(20) = str2double(get(param_20, 'String'));



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
   
    %# callback function
    function button_callback(src,ev)
        SWC = Data.SWC;                %Nt * Nz
        SoilT = Data.SoilT;            %Nt * Nz
        switch get(get(hBtnGrp,'SelectedObject'),'Tag')
            case 'Model_1',
                tic
                [out info] = DETECTmodel(Nt_run,SWC,SoilT,SWCant_Rm,SWCant_Rr,SoilTant,Gness,AtmPress,params_swp,Xrel,ic,d13Cr,d13Cm,dx,BC,params,params_lag);
                DayN=[(Nt_beg-1):0.25:Nt_end]';
                if (writefile==1)
                    save([pwd '/DETECT_Versions/Original_DETECT/Outputs/output_6hourly' num2str(ant) '.mat'],'out')
                    save([pwd '/DETECT_Versions/Original_DETECT/Outputs/info_6hourly' num2str(ant) '.mat'],'info')
                    myfile4 = fopen([pwd '/DETECT_Versions/Original_DETECT/Outputs/outCO2_6hourly' num2str(ant) '.csv'] ,'w+');
                    fprintf(myfile4,['DayN,' 'd0,' 'd1,' 'd2,' 'd3,' 'd4,' 'd5,' 'd6,' 'd7,' 'd8,' 'd9,' 'd10,' 'd11,' 'd12,' 'd13,' 'd14,' 'd15,' 'd16,' 'd17,' 'd18,' 'd19,' 'd20,' 'd21,' 'd22,' 'd23,' 'd24,' 'd25,' 'd26,' 'd27,' 'd28,' 'd29,' 'd30,' 'd31,' 'd32,' 'd33,' 'd34,' 'd35,' 'd36,' 'd37,' 'd38,' 'd39,' 'd40,'  'd41,' 'd42,' 'd43,' 'd44,' 'd45,' 'd46,' 'd47,' 'd48,' 'd49,' 'd50,'  'd51,' 'd52,' 'd53,' 'd54,' 'd55,' 'd56,' 'd57,' 'd58,' 'd59,' 'd60,'  'd61,' 'd62,' 'd63,' 'd64,' 'd65,' 'd66,' 'd67,' 'd68,' 'd69,' 'd70,'  'd71,' 'd72,' 'd73,' 'd74,' 'd75,' 'd76,' 'd77,' 'd78,' 'd79,' 'd80,'  'd81,' 'd82,' 'd83,' 'd84,' 'd85,' 'd86,' 'd87,' 'd88,' 'd89,' 'd90,'  'd91,' 'd92,' 'd93,' 'd94,' 'd95,' 'd96,' 'd97,' 'd98,' 'd99,' 'd100,\n']);
                    fclose(myfile4);
                    dlmwrite([pwd '/DETECT_Versions/Original_DETECT/Outputs/outCO2_6hourly' num2str(ant) '.csv'],[DayN out.CO2], '-append','precision', '%.3f');
                end
                res = ['Status: Done DETECT Original'];
                set(hEdit3, 'String',res)
                set(hEdit4, 'String',toc)
                set(hEdit4, 'Visible','on')
                set(hEdit5, 'Visible','on')
                if play_sound.Value == 1,
                    [y, Fs] = audioread([pwd '/Misc/finished.mp3']);
                    sound(y, Fs);
                end
            case 'Model_2',
                tic
                [out info] = Efficient_DETECTmodel(Nt_run,SWC,SoilT,SWCant_Rm,SWCant_Rr,SoilTant,Gness,AtmPress,params_swp,Xrel,ic,d13Cr,d13Cm,dx,BC,params,params_lag);
                DayN=[(Nt_beg-1):0.25:Nt_end]';
                if (writefile==1)
                    save([pwd '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/output_6hourly' num2str(ant) '.mat'],'out')
                    save([pwd '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/info_6hourly' num2str(ant) '.mat'],'info')
                    myfile4 = fopen([pwd '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/outCO2_6hourly' num2str(ant) '.csv'] ,'w+');
                    fprintf(myfile4,['DayN,' 'd0,' 'd1,' 'd2,' 'd3,' 'd4,' 'd5,' 'd6,' 'd7,' 'd8,' 'd9,' 'd10,' 'd11,' 'd12,' 'd13,' 'd14,' 'd15,' 'd16,' 'd17,' 'd18,' 'd19,' 'd20,' 'd21,' 'd22,' 'd23,' 'd24,' 'd25,' 'd26,' 'd27,' 'd28,' 'd29,' 'd30,' 'd31,' 'd32,' 'd33,' 'd34,' 'd35,' 'd36,' 'd37,' 'd38,' 'd39,' 'd40,'  'd41,' 'd42,' 'd43,' 'd44,' 'd45,' 'd46,' 'd47,' 'd48,' 'd49,' 'd50,'  'd51,' 'd52,' 'd53,' 'd54,' 'd55,' 'd56,' 'd57,' 'd58,' 'd59,' 'd60,'  'd61,' 'd62,' 'd63,' 'd64,' 'd65,' 'd66,' 'd67,' 'd68,' 'd69,' 'd70,'  'd71,' 'd72,' 'd73,' 'd74,' 'd75,' 'd76,' 'd77,' 'd78,' 'd79,' 'd80,'  'd81,' 'd82,' 'd83,' 'd84,' 'd85,' 'd86,' 'd87,' 'd88,' 'd89,' 'd90,'  'd91,' 'd92,' 'd93,' 'd94,' 'd95,' 'd96,' 'd97,' 'd98,' 'd99,' 'd100,\n']);
                    fclose(myfile4);
                    dlmwrite([pwd '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/outCO2_6hourly' num2str(ant) '.csv'],[DayN out.CO2], '-append','precision', '%.3f');
                end
                res = ['Status: Done DETECT Efficient'];
                set(hEdit3, 'String',res)
                set(hEdit4, 'String',toc)
                set(hEdit4, 'Visible','on')
                set(hEdit5, 'Visible','on')
                if play_sound.Value == 1,
                    [y, Fs] = audioread([pwd '/Misc/finished.mp3']);
                    sound(y, Fs);
                end
            case 'Model_3',
%                 SWC(:,:) = mean(SWC(:));
%                 SoilT(:,:) = mean(SoilT(:));
                tic
                [out info] = FAKE_SS_DETECTmodel(Nt_run,SWC,SoilT,SWCant_Rm,SWCant_Rr,SoilTant,Gness,AtmPress,params_swp,Xrel,ic,d13Cr,d13Cm,dx,BC,params,params_lag);
                DayN=[(Nt_beg-1):0.25:Nt_end]';
                if (writefile==1)
                    save([pwd '/DETECT_Versions/FAKE_SS_DETECT/FAKE_SS_Outputs/output_6hourly' num2str(ant) '.mat'],'out')
                    save([pwd '/DETECT_Versions/FAKE_SS_DETECT/FAKE_SS_Outputs/info_6hourly' num2str(ant) '.mat'],'info')
                    myfile4 = fopen([pwd '/DETECT_Versions/FAKE_SS_DETECT/FAKE_SS_Outputs/outCO2_6hourly' num2str(ant) '.csv'] ,'w+');
                    fprintf(myfile4,['DayN,' 'd0,' 'd1,' 'd2,' 'd3,' 'd4,' 'd5,' 'd6,' 'd7,' 'd8,' 'd9,' 'd10,' 'd11,' 'd12,' 'd13,' 'd14,' 'd15,' 'd16,' 'd17,' 'd18,' 'd19,' 'd20,' 'd21,' 'd22,' 'd23,' 'd24,' 'd25,' 'd26,' 'd27,' 'd28,' 'd29,' 'd30,' 'd31,' 'd32,' 'd33,' 'd34,' 'd35,' 'd36,' 'd37,' 'd38,' 'd39,' 'd40,'  'd41,' 'd42,' 'd43,' 'd44,' 'd45,' 'd46,' 'd47,' 'd48,' 'd49,' 'd50,'  'd51,' 'd52,' 'd53,' 'd54,' 'd55,' 'd56,' 'd57,' 'd58,' 'd59,' 'd60,'  'd61,' 'd62,' 'd63,' 'd64,' 'd65,' 'd66,' 'd67,' 'd68,' 'd69,' 'd70,'  'd71,' 'd72,' 'd73,' 'd74,' 'd75,' 'd76,' 'd77,' 'd78,' 'd79,' 'd80,'  'd81,' 'd82,' 'd83,' 'd84,' 'd85,' 'd86,' 'd87,' 'd88,' 'd89,' 'd90,'  'd91,' 'd92,' 'd93,' 'd94,' 'd95,' 'd96,' 'd97,' 'd98,' 'd99,' 'd100,\n']);
                    fclose(myfile4);
                    dlmwrite([pwd '/DETECT_Versions/FAKE_SS_DETECT/FAKE_SS_Outputs/outCO2_6hourly' num2str(ant) '.csv'],[DayN out.CO2], '-append','precision', '%.3f');
                end
                res = ['Status: Done DETECT FAKE SS'];
                set(hEdit3, 'String',res)
                set(hEdit4, 'String',toc)
                set(hEdit4, 'Visible','on')
                set(hEdit5, 'Visible','on')
%                 get(get(play_sound),'Value')
% play_sound.Value
                if play_sound.Value == 1,
                    [y, Fs] = audioread([pwd '/Misc/finished.mp3']);
                    sound(y, Fs);
                end
            case 'Model_4',
%                 SWC(:,:) = mean(SWC(:));
%                 SoilT(:,:) = mean(SoilT(:));
                tic
                [out info] = SS_DETECTmodel(Nt_run,SWC,SoilT,SWCant_Rm,SWCant_Rr,SoilTant,Gness,AtmPress,params_swp,Xrel,ic,d13Cr,d13Cm,dx,BC,params,params_lag);
                DayN=[(Nt_beg-1):0.25:Nt_end]';
                if (writefile==1)
                    save([pwd '/DETECT_Versions/SS_DETECT/SS_Outputs/output_6hourly' num2str(ant) '.mat'],'out')
                    save([pwd '/DETECT_Versions/SS_DETECT/SS_Outputs/info_6hourly' num2str(ant) '.mat'],'info')
                    myfile4 = fopen([pwd '/DETECT_Versions/SS_DETECT/SS_Outputs/outCO2_6hourly' num2str(ant) '.csv'] ,'w+');
                    fprintf(myfile4,['DayN,' 'd0,' 'd1,' 'd2,' 'd3,' 'd4,' 'd5,' 'd6,' 'd7,' 'd8,' 'd9,' 'd10,' 'd11,' 'd12,' 'd13,' 'd14,' 'd15,' 'd16,' 'd17,' 'd18,' 'd19,' 'd20,' 'd21,' 'd22,' 'd23,' 'd24,' 'd25,' 'd26,' 'd27,' 'd28,' 'd29,' 'd30,' 'd31,' 'd32,' 'd33,' 'd34,' 'd35,' 'd36,' 'd37,' 'd38,' 'd39,' 'd40,'  'd41,' 'd42,' 'd43,' 'd44,' 'd45,' 'd46,' 'd47,' 'd48,' 'd49,' 'd50,'  'd51,' 'd52,' 'd53,' 'd54,' 'd55,' 'd56,' 'd57,' 'd58,' 'd59,' 'd60,'  'd61,' 'd62,' 'd63,' 'd64,' 'd65,' 'd66,' 'd67,' 'd68,' 'd69,' 'd70,'  'd71,' 'd72,' 'd73,' 'd74,' 'd75,' 'd76,' 'd77,' 'd78,' 'd79,' 'd80,'  'd81,' 'd82,' 'd83,' 'd84,' 'd85,' 'd86,' 'd87,' 'd88,' 'd89,' 'd90,'  'd91,' 'd92,' 'd93,' 'd94,' 'd95,' 'd96,' 'd97,' 'd98,' 'd99,' 'd100,\n']);
                    fclose(myfile4);
                    dlmwrite([pwd '/DETECT_Versions/SS_DETECT/SS_Outputs/outCO2_6hourly' num2str(ant) '.csv'],[DayN out.CO2], '-append','precision', '%.3f');
                end
                res = ['Status: Done DETECT SS'];
                set(hEdit3, 'String',res)
                set(hEdit4, 'String',toc)
                set(hEdit4, 'Visible','on')
                set(hEdit5, 'Visible','on')
                if play_sound.Value == 1,
                    [y, Fs] = audioread([pwd '/Misc/finished.mp3']);
                    sound(y, Fs);
                end
             case 'Model_5',
                tic
                Rpath = 'D:\Programs\R\R_350\bin';
                RscriptFileName = 'D:\Research\DETECT_MatLab\DETECT_Versions\R_DETECT_NSS\DETECT_R_SCRIPT.R';
                RunRcode(RscriptFileName, Rpath); 

                res = ['Status: Done Rscript DETECT Efficient NSS'];
                set(hEdit3, 'String',res)
                set(hEdit4, 'String',toc)
                set(hEdit4, 'Visible','on')
                set(hEdit5, 'Visible','on')
                if play_sound.Value == 1,
                    [y, Fs] = audioread([pwd '/Misc/finished.mp3']);
                    sound(y, Fs);
                end
            otherwise, res = '';
        end
    end
            function button_callback2(hObject, eventdata, handles)
%                 Ns=[];
%                 DayN_RecoObs=[];
%                 RecoObs=[];
%                 DayN_RmObs=[];
%                 RmObs=[];
%                 DayN1 =[];
%                 DayN2 =[];
%                 precip_IN = [];
%                 precip=[];
%                 f1 =[];
%                 % hObject    handle to pushbutton1 (see GCBO)
%                 % eventdata  reserved - to be defined in a future version of MATLAB
%                 % handles    structure with handles and user data (see GUIDATA)
%                 display('Goodbye');
%                 h = figure;
%                 ax = axes('Parent',h, 'Units','normalized','Position', [0.1 0.1 0.8 0.6]);
% %                 plot(ax,magic(4))
%                  FAKE_SS_Rsoil_timeseries()
%                 
% %                 plot(0.5, 0.5, 'b*', 'MarkerSize', 15, 'LineWidth', 3);  % Adjust as needed.
% % grid on;
                    % clear all
                    % close all

                    Ns=1; 
                    DayN_RecoObs=[789 860 871 885 888 899 920 936 955 979] - (365*2);
                    RecoObs=[0.37 0.96 1.77 5.65 5.16 6.35 2.49 1.85 7.03 3.89]; %units: gC/m2/day
                    DayN_RmObs=[170 184 198 213 226 235 253 269];
                    RmObs=[1.48 1.22 0.52 0.46 2.54 1.61 1.31 1.01]; %units: gC/m2/day
                    DayN1=[92:274]';  %DoY for 2008/2009
                    DayN2=[457:639]'; %2008
                    precip_IN=csvread([pwd '/Plotting/Data_for_plots/Precip_PHACE_fromKevin_2007to2013_11Mar2015.csv'],1,1);
                    precip=repmat(precip_IN(DayN2,3),1,4);
                    f1=figure('Position',[10 40 1010 650]);           
                    for s=1:Ns
                        ah(s)=subplot(2,1,s)
                        if (s==1)
                            ant=['']; titlename=['(a) DETECT model without antecedent parameterization (\it{Ctrl} \rm{scenario)}'];
                        elseif (s==2)
                            ant=['_ant']; titlename=['(b) DETECT model with antecedent parameterization (\it{Ctrl-ant} \rm{scenario)}'];
                        end

                        switch get(get(hBtnGrp,'SelectedObject'),'Tag')
                        case 'Model_1',
                           load([pwd '/DETECT_Versions/Original_DETECT/Outputs/output_6hourly'  '.mat'],'out')  
                        case 'Model_2',
                            load([pwd '/DETECT_Versions/Efficient_DETECT/Efficient_Outputs/output_6hourly'  '.mat'],'out')
                        case 'Model_3', 
                            load([pwd '/DETECT_Versions/FAKE_SS_DETECT/FAKE_SS_Outputs/output_6hourly'  '.mat'],'out')
                        case 'Model_4',  
                            load([pwd '/DETECT_Versions/SS_DETECT/SS_Outputs/output_6hourly'  '.mat'],'out')
                        case 'Model_5',
                            % HERE NEED TO SET FOR R OUTPUT
                            load([pwd '/DETECT_Versions/R_DETECT_NSS/Outputs/output_6hourly'  '.mat'],'out')
                                otherwise, res = '';
                        end

                    %    load([pwd '/DETECT_Versions/FAKE_SS_DETECT/FAKE_SS_Outputs/output_6hourly'  '.mat'])
                       R_IN=out.CO2flux';
                       Rr_IN=sum(out.CO2fluxtype(1:2,:));
                       Rm_IN=sum(out.CO2fluxtype(3:4,:));
                       R_DETECT=mean(reshape(R_IN(1,1:732),4,183),1);
                       Rr_DETECT=mean(reshape(Rr_IN(1,1:732),4,183),1);
                       Rm_DETECT=mean(reshape(Rm_IN(1,1:732),4,183),1);     

%                         y = 5*sin(x);
%                         plot(ax,x,y)
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

                    % if (Ns==1)
                    %     saveas(f1,['Plots/FAKE_SS_Rsoil_timeseries_PHACEsite.jpg'])
                    % else
                    %     saveas(f1,['Rsoil_timeseries_withAnt_PHACEsite.jpg'])
                    % end




            end
end