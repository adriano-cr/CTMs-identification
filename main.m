clearvars
close all
clc

delete(gcp('nocreate'))
parpool('threads')

addpath(strcat(pwd,'\fnc\CTM'));
addpath(strcat(pwd,'\fnc\data_reader'));

%% Identification of CTM-s parameters
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


%% Choose desired code segments (true to execute, false to skip):
reader = true;         % Read and extract data from file
ctm = false;            % Identification of CTM parameters (no station)
ctm_s = true;           % Identification of CTM-s parameters (with station)
output_ctm = false;     % Output CTM data to CTM_param_out.xls
output_ctms = false;     % Output station data to CTM_param_out.xls

if(reader)
    %% 1. Data extraction
    % Extract the traffic data from the file obtained from Dexter NDW and
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

    station_params = CTMs_parameters(opt_CTMs);
    disp('CTM-s parameters identification done!')
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

    % Speed thresholds below which the traffic is assumed to be congested.
    % These are to be tuned for the problem at hand by looking at the plots
    % obtained in 'csv_DATEX_reader_v4'.

    opt_CTM.speed_th = [75 75 85 85 85 83 83 ...
        72 73 75 75 85 80 80 80];

    % The threshold used in the quantile regression, to be tuned for the
    % specific data in use. The vector must have as many elements as the
    % number of cells in the CTM-s model.

    opt_CTM.coeff_quantile = [0.98 0.5 0.99 0.98 0.75 0.75 0.75...
        0.80 0.80 0.75 0.75 0.75 0.75 0.75 0.75];

    [CTM_param,phi_1,last_phi] = CTM_identification(opt_CTM);
    disp('CTM parameters identification done!')
    disp('==============================')

    %% write output data
    path=strcat(pwd,'\fnc\extracted_data\CTM_param_out_A20_3lanes.xls');
    if(output_ctm)
        disp('==============================')
        fprintf('Saving CTM information in %s ...\n',path)

        ID = linspace(1,CTM_param.N,CTM_param.N).';
        L = round(CTM_param.len, 2);
        v = round(CTM_param.v_bar);
        w = round(CTM_param.w);
        q_max = round(CTM_param.q_max);
        rho_max = round(CTM_param.rho_max);
        t = CTM_param.T';
        tab = table(ID,L,v,w,q_max,rho_max,t);
        writetable(tab, path, 'Sheet','Cells parameters');
        writematrix(phi_1', path, 'Sheet','First Demand Real');
        writematrix(last_phi', path, 'Sheet','Last Demand Real');

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

    if(output_ctms)
        disp('==============================')
        fprintf('Saving station information in %s ...\n',path)

        beta = station_params.beta_avg;
        delta = station_params.delta_avg;
        occ = station_params.occupancy;
        s_s = station_params.flow_in;
        r_s = station_params.flow_out;

        tab = table(beta, delta);
        writetable(tab, path, 'Sheet','Station parameters');

        writetable(table(occ), path, 'Sheet','Occupancy station', 'WriteVariableNames', false);
        writetable(table(s_s), path, 'Sheet','Flow in station', 'WriteVariableNames', false);
        writetable(table(r_s), path, 'Sheet','Flow out station', 'WriteVariableNames', false);
    end
end

disp('... Finish!')
disp('==============================')