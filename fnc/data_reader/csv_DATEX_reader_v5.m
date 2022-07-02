function  out_structure = csv_DATEX_reader_v5(input_str,output_str,opt)
%% csv_DATEX_reader :
% The funciton that reads the raw traffic data file and then creates a
% structure as output with all the data nicely organized and ready to be
% elaborated. It prints also some rough plots of the data if required
% (the plot in this case are not automatically saved).
%
% INPUT :
%       - input_str : the string that define the input file where the raw
%       data are extracted
%       - output_str : the name of the file in which we want to save the
%       strucute that we compute
%       - extra.min_freq [1/min]: minimum frequency of the data
%       - extra.sensor_id : the name or code of the sensors
%       - opt.verbatim {0,1} : display or not information about the steps
%       performe
%       - opt.display {0,1} : plot or not the foundamental graphs associated
%       to the traffic
% OUTPUT :
%       - out_structure : the final structure with all the data
%                   The data structure is also saved in
%                   'fnc\data_reader\extracted_data\'output_str'.mat'

path1 = '\traffic_data';
addpath(genpath([pwd,path1]))
try
    %% Check number of inputs
    min_freq = opt.min_freq;
    % check if the frequency id an integer or not
    if ~mod(min_freq,1) == 0
        error('ERROR : "min_freq" has to be an integer')
    end

    %% Load data
    % Data obtained with th "volle" (full) structure
    filename = [input_str,'.csv'];

    disp('==============================')
    fprintf('Reading file: %s \n',filename)
    disp('==============================')

    % import the data file

    T = readtable(filename);

    %% Find sensor names

    [sensors_unique, sensors_index] = unique(T.naam_meetlocatie_mst);

    for i = 1:length(sensors_unique)
        sensors_unique(i) = erase(extractAfter(sensors_unique(i), 8), 'ra');
    end
    
    %% Extract useful data
    % extract from the whole dataset only the columns that interest us

    fprintf('Extracting useful data \n')
    T = (T(:,{'naam_meetlocatie_mst', 'start_meetperiode', 'eind_meetperiode', ...
        'gebruikte_minuten_intensiteit', 'gem_intensiteit', 'gem_snelheid', ...
        'totaal_aantal_rijstroken', 'rijstrook_rijbaan', 'start_locatie_latitude', ...
        'start_locatie_longitude', 'afstand_tot_vild_locatie'}));
    
    T = renamevars(T, {'naam_meetlocatie_mst', 'start_meetperiode', 'eind_meetperiode', ...
        'gebruikte_minuten_intensiteit', 'gem_intensiteit', 'gem_snelheid', ...
        'totaal_aantal_rijstroken', 'rijstrook_rijbaan', 'start_locatie_latitude', ...
        'start_locatie_longitude', 'afstand_tot_vild_locatie'}, ...
        {'sensor_id', 'start_time', 'end_time', ...
        'sample_time', 'avg_intensity', 'avg_speed', ...
        'num_lanes', 'lane_id', 'latitude', ...
        'longitude', 'position'});

    % collect all the data that are measured in the same time interval, they
    % are assumed consecutive

    fprintf('Reshaping the data \n')


    % create a temporary structure
    %sensor_sum(length(sensors_id)) = struct(); %preallocate space for speed-up

    T_extra = T;

    for i = 1:height(T)
        if cell2mat(extractAfter(T.lane_id(i), 4)) > 3
            T(i, :) = [];
        else
            T_extra(i, :) = [];
        end
    end
%             k = 1;
%          
%         for j = 1:3:length(sensor(i).veh_avg_speed)
%             veh1 = sensor(i).veh_number(j);
%             veh2 = sensor(i).veh_number(j+1);
%             veh3 = sensor(i).veh_number(j+2);
% 
%             total_veh = veh1 + veh2 + veh3;
% 
%             vel1 = sensor(i).veh_avg_speed(j);
%             vel2 = sensor(i).veh_avg_speed(j+1);
%             vel3 = sensor(i).veh_avg_speed(j+2);
% 
%             if total_veh == 0
%                 w_avg_speed = 0;
%             else
%                 w_avg_speed = (vel1 * veh1 + vel2 * veh2 + vel3 * veh3)/total_veh;
%             end
% 
%             sensor_sum(i).vehicle_number(k) = total_veh;
%             sensor_sum(i).vehicle_speed(k) = w_avg_speed;
% 
%             sensor_sum(i).ending_time(k) = sensor(i).ending_s_time(j);
%             sensor_sum(i).starting_time(k) = sensor(i).starting_s_time(j);
%             sensor_sum(i).sample_time(k) = sensor(i).time_sample(j);
%             sensor_sum(i).position(k) = sensor(i).location(j);
% 
%             k = k+1;
%         end
%     end

    sensor = sensor_sum;

    disp('==============================')


    %% Interpolate the data
    % if the minimum frequency is higher than the one of the
    % data we interpolate the data to attain the desired one
    if ~isempty(min_freq) && sensor(1).sample_time(1) > 1/min_freq

        for k = 1: length(sensors_id)
            xx = linspace(1, length(sensor(k).vehicle_number), length(sensor(k).vehicle_number));

            % every element in xx is made into "min_freq" many in yy
            yy = linspace(1,length(sensor(k).vehicle_number),length(sensor(k).vehicle_number)*min_freq);
            number_v = sensor(k).vehicle_number;
            interpolated_number_vv = interp1(xx,number_v,yy);
            % This is  the correct one because "interpolated_number_vv"
            % is already in [veh/h] ( changed wrt v1.0 )
            interpolated_number_vv = round(interpolated_number_vv);
            % interpolate the velocity
            speed_v = sensor(k).vehicle_speed;
            interpolated_speed_vv = interp1(xx,speed_v,yy);
            % assign the new values
            sensor(k).vehicle_speed = interpolated_speed_vv;
            % extend the other fields in "sensor"
            sensor(k).starting_time = repelem(sensor(k).starting_time,1,min_freq);
            sensor(k).ending_time = repelem(sensor(k).ending_time,1,min_freq);
            sensor(k).position = repelem(sensor(k).position,1,min_freq);

            % sample time in [h], from the site we have the data in
            % minutes hence we have to multiply 1/60 to achieve
            sample_interpolazione=sensor(k).sample_time(1)/min_freq;
            sensor(k).sample_time = sample_interpolazione/60*ones(1,length(yy));
            % Since the interpolated_number_vv is the number wrt to
            % hours, then we have to scale it wrt to the sample time.
            sensor(k).vehicle_number = interpolated_number_vv;
            %sensor(k).vehicle_number = interpolated_number_vv./60;
            %sensor(k).vehicle_number = interpolated_number_vv.*sensor(k).sample_time;
        end
    end

    %% Compute the  fundamental diagram
    % we already have the flow, we just need the density
    for j = 1:length(sensors_id)
        %flow = sensor(j).vehicle_number./sensor(j).sample_time; % [veh/h]
        %flow = sensor(j).vehicle_number.*60;
        flow=sensor(j).vehicle_number;
        density = flow./sensor(j).vehicle_speed;

        sensor(j).flow = flow;
        sensor(j).density = density;
    end
    % assign the output
    out_structure = sensor;
    %% Save the file
    save_file = [opt.path, output_str,'.mat'];
    save(save_file,'sensor')
    fprintf('6) Save the data in %s\n',save_file)
    disp('==============================')
    %% Plot

    if opt.display && length(sensor)>=4
        traffic_data_plot(sensor)
    end
catch ME
    keyboard
    rethrow(ME)
end
end

function traffic_data_plot(sensor)
%% Plot some data
last_fig_num = get(gcf,'Number');
n_row = 3; N = size(sensor,2);
for n = 1 : N
    % % % % % % % %
    figure(last_fig_num+1)
    subplot(n_row,ceil((N)/n_row),n)
    bar(sensor(n).vehicle_number)
    title_str1 = ['# vehicles (10s) (Sens. ',num2str(n),')'];
    title(title_str1)
    grid on

    % % % % % % % %
    figure(last_fig_num+2)
    subplot(n_row,ceil((N)/n_row),n)
    bar(sensor(n).vehicle_speed)
    grid on
    ax = gca;
    title_str2 = ['Avg. speed (Sens.',num2str(n),')'];
    title(title_str2)
    ax.YLim = [0,150];

    % % % % % % % %
    figure(last_fig_num+3)
    subplot(n_row,ceil((N)/n_row),n)
    scatter(sensor(n).density,sensor(n).flow)
    grid on
    title_str3 = ['Fundamental diagram (Sens. ',num2str(n),')'];
    title(title_str3)
end
end