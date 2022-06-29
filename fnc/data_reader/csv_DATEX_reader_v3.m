function  out_structure = csv_DATEX_reader_v3(input_str,output_str,opt,extra)
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

    min_freq = extra.min_freq;
    % check if the frequency id an integer or not
    if ~mod(min_freq,1) == 0
        error('ERROR : "min_freq" has to be an integer')
    end

    %% Check if the user wants to load other data
    % load a previously extracted data structure and immediately send it as
    % an output of the function
    if opt.load
        %         file_to_load = ['fnc\data_reader\extracted_data\',...
        %                         output_str,'.mat'];
        file_to_load = [output_str,'.mat'];
        load(file_to_load,'sensor')
        % assign the output
        out_structure = sensor;
        filename = [input_str,'.csv'];
        if opt.verbatim
            disp('==============================')
            fprintf('1) Loaded data in : %s \n',filename)
            disp('==============================')
        end
    else
        %% Load data
        % Data obtained with th "volle" (full) structure
        filename = [input_str,'.csv'];
        if opt.verbatim
            disp('==============================')
            fprintf('1) Use data in : %s \n',filename)
            disp('==============================')
        end
        % import in a cell
        import_raw = importdata(filename);
        cell_raw = import_raw.textdata;
        cell_raw(2:end,7) = num2cell(import_raw.data(:,1)); % fare dinamico
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
        if opt.verbatim
            fprintf('2) Construsction \n')
        end
        % complete all the fields with every row
        %
        data_field_name = fieldnames(data);
        for i = 1:numel(data_field_name)
            % get the data field
            field_i_name = char(data_field_name(i));
            % select the data in the row
            data.(field_i_name) = [cell_raw(2:end, i)];
        end
        if opt.verbatim
            fprintf('3) Created struct : %s \n','data')
        end
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
        if opt.verbatim
            fprintf('4) Organizing in struct sensor : %s \n', 'sensor')
        end
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
            sensor(j).location = data.afstand_tot_vild_locatie(sensor_index);
            
            % check corrupted data, hence zeros that should not be there
            % and due to sensor failure
            for k = 1 : length(sensor(j).veh_number)
                if(sensor(j).veh_number(k) == 0 )
                    if (k==1)
                        prev_num = 0;
                        prev_speed = 0;
                    else
                        prev_num = sensor(j).veh_number(k-1);
                        prev_speed = sensor(j).veh_avg_speed(k-1);
                    end
                    if (k==length(sensor(j).veh_number))
                        next_veh = 0;
                        next_speed = 0;
                    else
                        next_veh = sensor(j).veh_number(k+1);
                        next_speed = sensor(j).veh_avg_speed(k+1);
                    end
                    sensor(j).veh_number(k) = round((prev_num+next_veh)/2);
                    sensor(j).veh_avg_speed(k) = (prev_speed+next_speed)/2;
                    if(sensor(j).veh_avg_speed(k) < 0) 
                        sensor(j).veh_avg_speed(k) = 0;
                    end

                end
            end
        end

        % collect all the data that are measured in the same time interval, they
        % are assumed consecutive
        if opt.verbatim
            fprintf('5) Reshaping the data \n')
        end

        % create some temporary structure  fields
        sensor1 = sensor;
        [sensor1.('ending_time')] = sensor1.('ending_s_time');
        sensor1 = rmfield(sensor1,'ending_s_time');
        [sensor1.('starting_time')] = sensor1.('starting_s_time');
        sensor1 = rmfield(sensor1,'starting_s_time');
        [sensor1.('vehicle_number')] = sensor1.('veh_number');
        sensor1 = rmfield(sensor1,'veh_number');
        [sensor1.('vehicle_speed')] = sensor1.('veh_avg_speed');
        sensor1 = rmfield(sensor1,'veh_avg_speed');
        [sensor1.('position')] = sensor1.('location');
        sensor1 = rmfield(sensor1,'location');
        [sensor1.('sample_time')] = sensor1.('time_sample');
        sensor1 = rmfield(sensor1,'time_sample'); 
        
        % dovrebbe raggruppare le 3 corsie a parità di timestamp
        new_ind = 1;
        for k = 1:length(sensors_id)
            time_ind = '';
            for jj = 1 : length(sensor(k).starting_s_time)
                if strcmp(time_ind,sensor(k).ending_s_time(jj))
                    % it means that it is the same time
                    % create the weighted avg speed
                    total_veh = sensor1(k).vehicle_number(new_ind)+sensor(k).veh_number(jj);
                    v1 = sensor1(k).vehicle_speed(new_ind);
                    v2 = sensor(k).veh_avg_speed(jj);
                    % add the number of vehicles
                    sensor1(k).vehicle_number(new_ind) = total_veh;
                    sensor1(k).vehicle_speed(new_ind) = v1*sensor1(k).vehicle_number(new_ind)/total_veh ...
                        +v2*sensor(k).veh_number(jj)/total_veh;
                else
                    % there sample period that starts
                    % save the starting and ending time
                    sensor1(k).ending_time(new_ind) = sensor(k).ending_s_time(jj);
                    sensor1(k).starting_time(new_ind) = sensor(k).starting_s_time(jj);
                    sensor1(k).vehicle_number(new_ind) = sensor(k).veh_number(jj);
                    sensor1(k).vehicle_speed(new_ind) = sensor(k).veh_avg_speed(jj);
                    sensor1(k).sample_time(new_ind) = sensor(k).time_sample(jj);
                    sensor1(k).position(new_ind) = sensor(k).location(jj);
                    new_ind = new_ind+1;
                end
                time_ind = sensor(k).ending_s_time(jj);
            end
            sensor1(k).ending_time(new_ind:end) = [];
            sensor1(k).starting_time(new_ind:end) = [];
            sensor1(k).vehicle_number(new_ind:end) = [];
            sensor1(k).vehicle_speed(new_ind:end) = [];
            sensor1(k).sample_time(new_ind:end) = [];
            sensor1(k).position(new_ind:end) = [];
            new_ind = 1;
            % make all column vectors
            sensor1(k).ending_time = sensor1(k).ending_time';
            sensor1(k).starting_time = sensor1(k).starting_time';
            sensor1(k).position = sensor1(k).position';
        end

        sensor = sensor1;
        if opt.verbatim
            disp('==============================')
        end

        %% Interpolate the data

        if ~isempty(min_freq) && sensor(1).sample_time(1) > 1/min_freq
            % if the minimum frequency is higher than the one of the
            % data we interpolate the data to attain the desired one
            for k = 1: length(sensors_id)
                xx = 1:length(sensor(k).vehicle_number);
                % every element in xx is made into "min_freq" many in yy
                yy = 1:1/min_freq:length(sensor(k).vehicle_number);
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
                % delete extra elements
                sensor(k).starting_time(end-2:end) = [];
                sensor(k).ending_time(end-2:end) = [];
                sensor(k).position(end-2:end) = [];

                % sample time in [h], from the site we have the data in
                % minutes hence we have to multiply 1/60 to achieve
                sensor(k).sample_time = sensor(k).sample_time(1)/min_freq*ones(1,length(yy))*(1/60);
                % Since the interpolated_number_vv is the number wrt to
                % hours, then we have to scale it wrt to the sample time.
                sensor(k).vehicle_number = interpolated_number_vv.*sensor(k).sample_time;
            end
        end

        %% Compute the  fundamental diagram
        % we already have the flow, we just need the density
        for j = 1:length(sensors_id)
            % flow = vehicle_num / sample time
            % number every hour ( changed wrt v1.0 )
            flow = sensor(j).vehicle_number./sensor(j).sample_time; % [veh/h] 
            density = flow./sensor(j).vehicle_speed;
            sensor(j).flow = flow;
            sensor(j).density = density;
        end
        % assign the output
        out_structure = sensor;
        %% Save the file
        %save_file = ['C:\A_Tesi\CTM-identification\fnc\data_reader\extracted_data', output_str,'.mat'];
        save_file = ['C:/Users/adria/Documents/Uni/LM II anno/Tesi/CTM-identification/fnc/data_reader/extracted_data/', output_str,'.mat'];
        %save_file = ['H:\Il mio Drive\Tesi magistrale\CTM-identification\fnc\data_reader\extracted_data', output_str,'.mat'];
        save(save_file,'sensor')
        if opt.verbatim
            fprintf('6) Save the data in %s\n',save_file)
            disp('==============================')
        end
    end
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
    title_str1 = ['# vehicles (Sens. ',num2str(n),')'];
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