clearvars
close all
clc

addpath(strcat(pwd,'\fnc\CTM'));
addpath(strcat(pwd,'\fnc\data_reader'));

reader = false;
ctm_s = true;
ctm = false;
output_data = false; % Output identified data to CTM_param_out.xls

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
% The data on traffic flows in the Netherlands are obtained from
% https://dexter.ndwcloud.nu/home
% The user is assumed to select the the sensors and then save a .csv file.
% This has to be saved in the folder '>fnc>data_reader>traffic_data'
% The details on the data used in this example can be found in detail in
% the papaer "Highway Traffic Control via Smart
% e-Mobility - Part II : Case Study"


if(reader)
    %% 1. Data extraction
    % Extract the data about traffic from the file 'intensiteit-snelheid-export'
    % stored in the folder: >fnc>data_reader>traffic_data

    % Plot option for the data obtained graphs (1 to turn on, 0 to turn off)
    opt_DATEX.display = 1;

    % The number of the lane used to enter/exit the service station.
    % e.g.: "lane3" for the rightmost lane of a 3 lanes road
    % Leave as "" if this functionality is not needed.
    opt_DATEX.laneSS = "";

    csv_DATEX_reader_v4('A2-southbound-station',opt_DATEX);
    disp('Reading done!')
    disp('==============================')
end

if(ctm_s)
    %% 2. CTMs param identification
    % Plot option for the data obtained graphs (1 to turn on, 0 to turn off)
    opt_CTMs.display = 1;

    % The ID of the sensor right before and after the service station of a
    % stretch, e.g.: "101", "201", etc.
    opt_CTMs.id_sensor_input = "522";
    opt_CTMs.id_sensor_output = "535";

    CTMs_parameters(opt_CTMs);
    disp('CTMs parameters identification done!')
    disp('==============================')
end

if(ctm)
    %% 3. CTM param identification
    % Identify the parameters of the CTM model
    % and the input flow 'phi_1'

    % Plot option for the identification related graphs (1 to turn on, 0 to turn off)
    opt_CTM.disp = 1;

    % Tolerance for the identification of outliers in the fundamental graph
    opt_CTM.tolerance = 500;

    % Speed threshold below which the vehicles are assumed to be into a congestion.
    % This is to be tuned for the problem at hand by looking at the plots
    % obtained in 'csv_DATEX_reader_v4'.

    opt_CTM.speed_th = [88 85 80 85 85 83 83 ...
        72 73 84 82 85 80 ];

    % The threshold used in the quantile regression.
    % These are to be tuned for the particular data used.
    % The vector has to be as long as the number of cells in the CTM-s model.

    opt_CTM.coeff_quantile = [0.98 0.98 0.75 0.75 0.75 0.75 0.75...
        0.80 0.80 0.75 0.75 0.75 0.75 ];

    [CTM_param,phi_1,last_phi] = CTM_identification(opt_CTM);
    disp('CTM parameters identification done!')
    disp('==============================')
    %% write output data
    if(output_data)
        path=strcat(pwd,'\fnc\extracted_data\CTM_param_out.xls');
        disp('==============================')
        fprintf('Saving information in %s ...\n',path)
        
        ID = linspace(1,CTM_param.N,CTM_param.N).';
        L = round(CTM_param.len, 2);
        v = round(CTM_param.v_bar);
        w = round(CTM_param.w);
        q_max = round(CTM_param.q_max);
        rho_max = round(CTM_param.rho_max);
        t = CTM_param.T';
        tabella = table(ID,L,v,w,q_max,rho_max,t);
        writetable(tabella, path, 'Sheet','Cells parameters');

        phi_1=round(phi_1*1.5);
        phi_1=table(phi_1.');

        last_phi=round(last_phi*1.5);
        last_phi=table(last_phi.');

        d = 0:seconds(10):hours(24);
        d = d(2:end);
        phi_smooth = table2timetable(phi_1, 'RowTimes',d');
        phi_smooth = smoothdata(phi_smooth,"sgolay","SmoothingFactor",0.15,"Degree",4);
        phi_smooth = timetable2table(phi_smooth);
        phi_smooth = phi_smooth(:,"Var1");

        phi_smooth_last = table2timetable(phi_1, 'RowTimes',d');
        phi_smooth_last = smoothdata(phi_smooth_last,"sgolay","SmoothingFactor",0.15,"Degree",4);
        phi_smooth_last = timetable2table(phi_smooth_last);
        phi_smooth_last = phi_smooth_last(:,"Var1");
  
        writetable(phi_smooth, path, 'Sheet','First Demand Smooth', 'WriteVariableNames', false);
        writetable(phi_smooth_last, path, 'Sheet','Last Demand Smooth', 'WriteVariableNames', false);

    end
end
disp('... Finish!')
disp('==============================')