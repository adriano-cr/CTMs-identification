clearvars
close all
clc
p = genpath('fnc');
addpath(p);

%% Identification of CTM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original code by: Carlo Cenedese (ccenedese@ethz.ch)
% Affiliation: ETH zurich - Automatic Control Laboratory
% Based on the theory developed in "Highway Traffic Control via Smart
% e-Mobility" submitted to IEEE - TRANSACTION ON INTELLIGENT TRANSPORTATION
% SYSTEMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forked and reworked for dynamical use with CTM-s models by:
% Davide Spalenza (davide.spalenza01@universitadipavia.it) &
% Adriano Cotta Ramusino (adriano.cottaramusino01@universitadipavia.it)
% Code version: v5.2
% Date: 08/09/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preliminaries
% The data on traffic flows in the netherlands are obtained from
% https://dexter.ndwcloud.nu/home
% The user is assumed to select the the sensors and then save a .csv file.
% This has to be saved in the folder '>fnc>data_reader>traffic_data'
% The details on the data used in this example can be found in details in
% the papaer "Highway Traffic Control via Smart
% e-Mobility - Part II : Case Study"


%% 1. Data extraction
% Extract the data about traffic from the file 'intensiteit-snelheid-export' stored in
% the folder: >fnc>data_reader>traffic_data

opt_DATEX.path = 'C:\A_Tesi\CTMs-identification\fnc\data_reader\extracted_data';
%opt_DATEX.path = 'C:/Users/adria/Documents/Uni/LM II anno/Tesi/CTMs-identification/fnc/data_reader/extracted_data/';
%opt_DATEX.path = 'H:\Il mio Drive\Tesi magistrale\CTMs-identification\fnc\data_reader\extracted_data';

% Plot option for the data obtained graphs (1 to turn on, 0 to turn off)
opt_DATEX.display = 0;

% The number of the lane used to enter/exit the service station.
% e.g.: "lane3" for the rightmost lane of a 3 lanes road
% Leave as "" if this functionality is not needed.
opt_DATEX.laneSS = "";

% The ID of the sensor right before and after the service station of a
% stretch, e.g.: "101", "201", etc.
opt_DATEX.id_sensor_input = "522";
opt_DATEX.id_sensor_output = "535";

data = csv_DATEX_reader_v4('A2-southbound-station','data_structure4_v2',opt_DATEX);

%% 2. CTM param identification
% Identify the parameters of the CTM model
% and the input flow 'phi_1'

% Plot option for the identification related graphs (1 to turn on, 0 to turn off)
opt_identification.disp = 1;

% Tolerance for the identification of outliers in the fundamental graph
opt_identification.tolerance = 500;

% Speed threshold below which the vehicles are assumed to be into a congestion.
% This is to be tuned for the problem at hand by looking at the plots
% obtained in 'csv_DATEX_reader_v4'.

opt_identification.speed_th = [88 85 80 85 85 83 83 ...
                               72 73 84 82 85 80 ];


% The threshold used in the quantile regression.
% These are to be tuned for the particular data used.
% The vector has to be as long as the number of cells in the CTM-s model.

opt_identification.coeff_quantile = [0.98 0.98 0.75 0.75 0.75 0.75 0.75...
                                     0.80 0.80 0.75 0.75 0.75 0.75 ];

% Output identified data to CTM_param_out.xls (1 to turn on, 0 to turn off)
output_data = 1;

[CTM_param,phi_1] = CTM_identification(data,opt_identification);

%% write output data
if(output_data > 0)
    fprintf('7) Saving information in CTM_param_out.xls... \n')
    ID = linspace(1,CTM_param.N,CTM_param.N).';
    L = round(CTM_param.len, 2);
    v = round(CTM_param.v_bar);
    w = round(CTM_param.w);
    q_max = round(CTM_param.q_max);
    rho_max = round(CTM_param.rho_max);
    T = table(ID,L,v,w,q_max,rho_max);
    writetable(T, "CTM_param_out.xls", 'Sheet',1);

    phi_1=round(phi_1*1.5);
    phi_1=table(phi_1.');

    d = 0:seconds(10):hours(24);
    d = d(2:end);
    phi_smooth = table2timetable(phi_1, 'RowTimes',d');
    phi_smooth = smoothdata(phi_smooth,"sgolay","SmoothingFactor",0.15,"Degree",4);
    phi_smooth = timetable2table(phi_smooth);
    phi_smooth = phi_smooth(:,"Var1");

    writetable(phi_1, "CTM_param_out.xls", 'Sheet',2, 'WriteVariableNames', false);
    writetable(phi_smooth, "CTM_param_out.xls", 'Sheet',3, 'WriteVariableNames', false);
end
fprintf('7) ... done!\n')

