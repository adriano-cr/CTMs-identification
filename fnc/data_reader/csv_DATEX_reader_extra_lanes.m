function  out_structure = csv_DATEX_reader_extra_lanes(input_str,output_str,opt)
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
    fprintf('1) Use data in : %s \n',filename)
    disp('==============================')

    % import in a cell
    import_raw = importdata(filename, ';');
    cell_raw = import_raw.textdata;
    cell_raw(2:end,12) = num2cell(import_raw.data(:,1)); % fare dinamico
    % create an empty structure
    data = struct();

    %% Header
    % the first row is the header,
    % thus it is used to create the structure field
    header = cell_raw(1,:);

    % create all the fields
    for j = 1 : length(header)
        field_name = char(header(1,j));
        if ~isnan(str2double(field_name(1)))
            % to avoid error due to first el a number
            %field_name = ['a', field_name];
        end
        data.(field_name) = [];
        % setfield(data,field_name,[]);
    end

    %% Fill the data fields

    fprintf('2) Construsction \n')

    % complete all the fields with every row
    %
    data_field_name = fieldnames(data);
    for i = 1:numel(data_field_name)
        % get the data field
        field_i_name = char(data_field_name(i));
        % select the data in the row
        data.(field_i_name) = [cell_raw(2:end, i)];
    end

    fprintf('3) Created struct : %s \n','data')

    sensors_raw = unique(data.naam_meetlocatie_mst);

    %% Find sensor names
    % sensors_id = char(1, length(sensors_raw));
    for i = 1:length(sensors_raw)
        sensors_id(i) = erase(extractAfter(string(sensors_raw(i)), 8), 'ra');
    end
    %sensors_id = fliplr(sensors_id);

    %% Extract useful data
    % extract from the whole data only the ones that interest us
    % find different inde associated to the different sensors

    fprintf('4) Organizing in struct sensor : %s \n', 'sensor')

    sensor(length(sensors_id)) = struct(); %preallocate space for speed-up
    for j = 1:length(sensors_id)

        check_sensor = strfind(data.naam_meetlocatie_mst,sensors_id(j));
        sensor_index = zeros(length(check_sensor),1);

        for i = 1:length(check_sensor)
            sensor_index(i) = ~isempty(cell2mat(check_sensor(i)));
        end
        sensor_index = logical(sensor_index);

        sensor(j).id = sensors_id(j);
        sensor(j).veh_number = str2double(data.gem_intensiteit(sensor_index)); % the intensity of vehicles over one h
        sensor(j).veh_avg_speed = str2double(data.gem_snelheid(sensor_index));
        sensor(j).time_sample = str2double(data.gebruikte_minuten_intensiteit(sensor_index)); %sample time is not always consistent even though the most cases it is
        sensor(j).ending_s_time = data.eind_meetperiode(sensor_index);
        sensor(j).starting_s_time = data.start_meetperiode(sensor_index);
        sensor(j).lane = data.rijstrook_rijbaan(sensor_index);
    end

 

    fprintf('5) Reshaping the data \n')
    
    % check corrupted data, hence zeros that should not be there
    % and due to sensor failure
    for j = 1:length(sensors_id)
        for k = 1 : length(sensor(j).veh_number)
                if(sensor(j).veh_avg_speed(k)<0)
                    sensor(j).veh_avg_speed(k) = 0;
                end
        end
    end

    % collect all the data that are measured in the same time interval, they
    % are assumed consecutive
    % create a temporary structure


    disp('==============================')


    %% Interpolate the data
    % if the minimum frequency is higher than the one of the
    % data we interpolate the data to attain the desired one
    if ~isempty(min_freq) && sensor(1).time_sample(1) > 1/min_freq

        for k = 1: length(sensors_id)
            xx = 1:length(sensor(k).veh_number);
            % every element in xx is made into "min_freq" many in yy
            yy = 1:1/min_freq:length(sensor(k).veh_number);
            number_v = sensor(k).veh_number;
            interpolated_number_vv = interp1(xx,number_v,yy);
            % This is  the correct one because "interpolated_number_vv"
            % is already in [veh/h] ( changed wrt v1.0 )
            interpolated_number_vv = round(interpolated_number_vv);
            % interpolate the velocity
            speed_v = sensor(k).veh_avg_speed;
            interpolated_speed_vv = interp1(xx,speed_v,yy);
            % assign the new values
            sensor(k).veh_avg_speed = interpolated_speed_vv;
            % extend the other fields in "sensor"
            sensor(k).starting_s_time = repelem(sensor(k).starting_s_time,1,min_freq);
            sensor(k).ending_s_time = repelem(sensor(k).ending_s_time,1,min_freq);
%             sensor(k).latitude = repelem(sensor(k).latitude, 1, min_freq);
%             sensor(k).longitude = repelem(sensor(k).longitude, 1, min_freq);
            sensor(k).lane = repelem(sensor(k).lane, 1, min_freq);

            % sample time in [h], from the site we have the data in
            % minutes hence we have to multiply 1/60 to achieve
            sensor(k).time_sample = sensor(k).time_sample(1)/min_freq*ones(1,length(yy))*(1/60);
            % Since the interpolated_number_vv is the number wrt to
            % hours, then we have to scale it wrt to the sample time.
            %sensor(k).vehicle_number = interpolated_number_vv;
            %sensor(k).vehicle_number = interpolated_number_vv./60;
            sensor(k).veh_number = interpolated_number_vv.*sensor(k).time_sample;
        end
    end

    %% Compute the  fundamental diagram
    % we already have the flow, we just need the density
    for j = 1:length(sensors_id)
        flow = sensor(j).veh_number./sensor(j).time_sample; % [veh/h]
 
        %flow = (sensor(j).vehicle_number./sensor(j).sample_time)./3; % [veh/h]
        density = flow./sensor(j).veh_avg_speed;

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
    bar(sensor(n).veh_number)
    title_str1 = ['# vehicles (10s)-sens',sensor(n).id];
    title(title_str1)
    grid on

    % % % % % % % %
    figure(last_fig_num+2)
    subplot(n_row,ceil((N)/n_row),n)
    bar(sensor(n).veh_avg_speed)
    grid on
    ax = gca;
    title_str2 = ['Avg. speed-sens',sensor(n).id];
    title(title_str2)
    ax.YLim = [0,150];

    % % % % % % % %
    figure(last_fig_num+3)
    subplot(n_row,ceil((N)/n_row),n)
    scatter(sensor(n).density,sensor(n).flow)
    grid on
    title_str3 = ['Fundamental diagram-sens',sensor(n).id];
    title(title_str3)
end
end