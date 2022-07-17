clearvars
close all
clc
p = genpath('fnc');
addpath(p);

%% Identification of CTM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Carlo Cenedese (ccenedese@ethz.ch)
% Affiliation: ETH zurich - Automatic Control Laboratory
% Code version: v2.0
% Date: 21/01/2022
% Based on the theory developed in "Highway Traffic Control via Smart
% e-Mobility" submitted to IEEE - TRANSACTION ON INTELLIGENT
% TRANSPORTATION SYSTEMS available at

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preliminaries
% The data on traffic flows in the netherlands are obtained from
% https://dexter.ndwcloud.nu/home
% The user is assumed to select the the sensors and then dave a .cvs data.
% This has to be saved in the folder '>fnc>data_reader>traffic_data'
% The details on the data used in this example can be found in details in
% the papaer "Highway Traffic Control via Smart
% e-Mobility - Part II : Case Study"


%% 1. Data extraction
% Extract the data about traffic from the file 'intensiteit-snelheid-export' stored in
% the folder: >fnc>data_reader>traffic_data

% plot the data obtained graphs
opt_DATEX.display = 1;

opt_DATEX.path = 'C:\A_Tesi\CTMs-identification\fnc\data_reader\extracted_data';
%opt_DATEX.path = 'C:/Users/adria/Documents/Uni/LM II anno/Tesi/CTMs-identification/fnc/data_reader/extracted_data/';
%opt_DATEX.path = 'H:\Il mio Drive\Tesi magistrale\CTMs-identification\fnc\data_reader\extracted_data';

% minimum frequency per minute at which you want the
% data. Performs an interpolation if the data are sampled at a lower rate
opt_DATEX.min_freq = 6;

opt_DATEX.laneSS = ""; 
opt_DATEX.id_sensor_input = "522";
opt_DATEX.id_sensor_output = "535";

data = csv_DATEX_reader_v4('A2-southbound-station','data_structure4_v2',opt_DATEX);

%% 2. CTM param identification
% Identify the parameters of the CTM model
% and the input flow 'phi_1'

% plot the figures realtive to the CTM identification, mainly the
% fundamental diagram
opt_identification.disp = 1;

opt_identification.tolerance = 500;

% threshold below which the vehicles are assumed to be into a congestion.
% This was tuned for the problem at hand by looking at the plots attained
% in 'cvs_DATEX_reader' regarding the velocity.

opt_identification.speed_th = [85 85 80 90 90 90 90 ...
                               80 85 90 90 85 80 85];


% The threshold used in the quantile regression, these are tuned for the
% particular date used in this example. The vector has to be as long as the
% number of cells in teh CTM model.

opt_identification.coeff_quantile = [0.82 0.75 0.80 0.90 0.90 0.90 0.95...
                                     0.80 0.90 0.90 0.90 0.90 0.90 0.90];
% 
% 
% [CTM_param,phi_1] = CTM_identification(data,opt_identification);
% 
% %% write output data
% ID=linspace(1,CTM_param.N,CTM_param.N).';
% L=CTM_param.len;
% v=CTM_param.v_bar;
% w=CTM_param.w;
% q_max=CTM_param.q_max;
% rho_max=CTM_param.rho_max;
% T = table(ID,L,v,w,q_max,rho_max);
% writetable(T, "CTM_param_out.xls", 'Sheet',1);
% 
% phi_1=table(phi_1.');
% writetable(phi_1, "CTM_param_out.xls", 'Sheet',2, 'WriteVariableNames', false);