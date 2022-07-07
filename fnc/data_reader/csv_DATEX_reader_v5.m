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
    fprintf('1) Use data in : %s \n',filename)
    disp('==============================')

    % import in a cell
    import_raw = importdata(filename, ';');
    cell_raw = import_raw.textdata;
    cell_raw(2:end,9) = num2cell(import_raw.data(:,1)); % fare dinamico
    cell_raw(2:end,10) = num2cell(import_raw.data(:,2)); % fare dinamico

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
    lane_ss = "lane6";
    sensor(length(sensors_id)) = struct(); %preallocate space for speed-up

    for j = 1:length(sensors_id)

        check_sensor = strfind(erase(extractAfter(string(data.naam_meetlocatie_mst), 8), 'ra'),sensors_id(j));
        sensor_index = zeros(length(check_sensor),1);

        for i = 1:length(check_sensor)
            sensor_index(i) = ~isempty(cell2mat(check_sensor(i)));
        end
        sensor_index = logical(sensor_index);
        sensor(j).lane = data.rijstrook_rijbaan(sensor_index);
        sensor(j).id = sensors_id(j);
        sensor(j).veh_number = str2double(data.gem_intensiteit(sensor_index)); % the intensity of vehicles over one h
        sensor(j).veh_avg_speed = str2double(data.gem_snelheid(sensor_index));
        sensor(j).time_sample = str2double(data.gebruikte_minuten_intensiteit(sensor_index)); %sample time is not always consistent even though the most cases it is
        sensor(j).ending_s_time = data.eind_meetperiode(sensor_index);
        sensor(j).starting_s_time = data.start_meetperiode(sensor_index);
        sensor(j).latitude = data.start_locatie_latitude(sensor_index);
        sensor(j).longitude = data.start_locatie_longitude(sensor_index);
    end

    kk=1;
    jj=1;
    flag=false;
    flag2=false;
    for j = 1:length(sensors_id)
        k=1;
        p=1;
        for i = 1 : length(sensor(j).lane)
            app = sensor(j).lane(i);
            if(strcmp(lane_ss, app))
                sensor_ss(kk).id = sensor(j).id;
                sensor_ss(kk).veh_number(k) = sensor(j).veh_number(i);
                sensor_ss(kk).veh_avg_speed(k) = sensor(j).veh_avg_speed(i);
                if(sensor_ss(kk).veh_avg_speed(k)<0)
                    sensor_ss(kk).veh_avg_speed(k)=0;
                end
                sensor_ss(kk).time_sample(k) = sensor(j).time_sample(i); %sample time is not always consistent even though the most cases it is
                sensor_ss(kk).ending_s_time(k) = sensor(j).ending_s_time(i);
                sensor_ss(kk).starting_s_time(k) = sensor(j).starting_s_time(i);
                sensor_ss(kk).latitude(k) =  sensor(j).latitude(i);
                sensor_ss(kk).longitude(k) = sensor(j).longitude(i);
                sensor_ss(kk).lane(k) =  sensor(j).lane(i);
                k=k+1;
                flag=true;
            else
                sensor_main(jj).id = sensor(j).id;
                sensor_main(jj).veh_number(p) = sensor(j).veh_number(i);
                sensor_main(jj).veh_avg_speed(p) = sensor(j).veh_avg_speed(i);
                if(sensor_main(jj).veh_avg_speed(p)<0)
                    sensor_main(jj).veh_avg_speed(p)=0;
                end
                sensor_main(jj).time_sample(p) = sensor(j).time_sample(i); %sample time is not always consistent even though the most cases it is
                sensor_main(jj).ending_s_time(p) = sensor(j).ending_s_time(i);
                sensor_main(jj).starting_s_time(p) = sensor(j).starting_s_time(i);
                sensor_main(jj).latitude(p) =  sensor(j).latitude(i);
                sensor_main(jj).longitude(p) = sensor(j).longitude(i);
                sensor_main(jj).lane(p) =  sensor(j).lane(i);
                p=p+1;
                flag2=true;
            end
        end
        if flag
            kk=kk+1;
            flag=false;
        end
        if flag2
            jj=jj+1;
            flag2=false;
        end
    end

    fprintf('5) Reshaping the data \n')
    % collect all the data that are measured in the same time interval, they
    % are assumed consecutive
    % create a temporary structure
    sensor_sum(length(sensors_id)) = struct(); %preallocate space for speed-up

    for i = 1:length(sensors_id)
        k = 1;
        for j = 1:5:length(sensor_main(i).starting_s_time)
            veh1 = sensor_main(i).veh_number(j);
            veh2 = sensor_main(i).veh_number(j+1);
            veh3 = sensor_main(i).veh_number(j+2);
            veh4 = sensor_main(i).veh_number(j+3);
            veh5 = sensor_main(i).veh_number(j+4);

            total_veh = veh1 + veh2 + veh3+veh4+veh5;

            vel1 = sensor_main(i).veh_avg_speed(j);
            vel2 = sensor_main(i).veh_avg_speed(j+1);
            vel3 = sensor_main(i).veh_avg_speed(j+2);
            vel4 = sensor_main(i).veh_avg_speed(j+3);
            vel5 = sensor_main(i).veh_avg_speed(j+4);
            if total_veh == 0
                w_avg_speed = 0;
            else
                w_avg_speed = (vel1 * veh1 + vel2 * veh2 + vel3 * veh3 + veh4*vel4 + veh5*vel5)/total_veh;
            end
            sensor_sum(i).vehicle_number(k) = total_veh;
            sensor_sum(i).vehicle_speed(k) = w_avg_speed;
            sensor_sum(i).latitude(k) = sensor_main(i).latitude(j);
            sensor_sum(i).longitude(k) = sensor_main(i).longitude(j);
            sensor_sum(i).lane(k) = sensor_main(i).lane(j);
            sensor_sum(i).id = sensor_main(i).id;
            sensor_sum(i).ending_time(k) = sensor_main(i).ending_s_time(j);
            sensor_sum(i).starting_time(k) = sensor_main(i).starting_s_time(j);
            sensor_sum(i).sample_time(k) = sensor_main(i).time_sample(j);
            k = k+1;
        end
    end

    %% "interpolate"

    
    %sensor_interp(int8(length(sensor_sum(1).starting_time)/6))=struct();
    sensor_interp(length(sensors_id))=struct();
    sensor_interp.sensor(1440)=struct();

   for i = 1:length(sensors_id)
    for j=1:length(sensor_sum(i).starting_time)

        hour = string(extractAfter(sensor_sum(i).starting_time(j), 11));


        check_hour = strfind(sensor_sum(i).starting_time,hour);
        hour_index = zeros(length(check_hour),1);

        for k = 1:length(check_hour)
            hour_index(k) = ~isempty(cell2mat(check_hour(k)));
        end

        hour_index = logical(hour_index);

        sensor_interp(i).starting_time = hour;
        sensor_interp(i).vehicle_number = sensor_sum(i).vehicle_number(hour_index);
        sensor_interp(i).vehicle_speed = sensor_sum(i).vehicle_speed(hour_index);
    end
   end

    %% Compute the  fundamental diagram
    % we already have the flow, we just need the density
    for j = 1:length(sensor_interp)
        flow = sensor_interp(j).vehicle_number; % [veh/h]

        density = flow./sensor_interp(j).vehicle_speed;
        sensor_interp(j).flow = flow;
        sensor_interp(j).density = density;
    end

    % assign the output
    out_structure = sensor_sum;
    %% Save the file
    save_file = [opt.path, output_str,'.mat'];
    save(save_file,'sensor')
    fprintf('6) Save the data in %s\n',save_file)
    disp('==============================')
    %% Plot
    traffic_data_plot(sensor_interp)

catch ME
    keyboard
    rethrow(ME)
end
end

function traffic_data_plot(sensor_sum)
%% Plot some data
last_fig_num = get(gcf,'Number');
n_row = 3; N = size(sensor_sum,2);
for n = 1 : N
    % % % % % % % %
    figure(last_fig_num+1)
    subplot(n_row,ceil((N)/n_row),n)
    bar(sensor_interp(n).vehicle_number)
    title_str1 = ['# vehicles (10s)-sens',sensor_sum(n).id];
    title(title_str1)
    grid on

    % % % % % % % %
    figure(last_fig_num+2)
    subplot(n_row,ceil((N)/n_row),n)
    bar(sensor_interp(n).vehicle_speed)
    grid on
    ax = gca;
    title_str2 = ['Avg. speed-sens',sensor_sum(n).id];
    title(title_str2)
    ax.YLim = [0,150];

    % % % % % % % %
    figure(last_fig_num+3)
    subplot(n_row,ceil((N)/n_row),n)
    scatter(sensor_interp(n).density,sensor_interp(n).flow)
    grid on
    title_str3 = ['Fundamental diagram-sens',sensor_sum(n).id];
    title(title_str3)
end
end