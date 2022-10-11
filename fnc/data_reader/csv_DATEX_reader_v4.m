function  [] = csv_DATEX_reader_v4(input_str,opt)
%% csv_DATEX_reader :
% Function that reads the raw traffic data file and then creates a
% structure as output with all the data nicely organized and ready to be
% elaborated. It also displays some plots of the data if required
% (the plots in this case are not automatically saved).
%
% INPUT :
%       - input_str: the string that defines the input file where the raw
%                    data are stored
%       - output_str: the name of the file in which the computed structure
%                    is to be saved
%       - opt.display {0,1}: 1 to plot the foundamental graphs associated
%                    with the traffic data, 0 to omit
%       - opt.laneSS: name of service station lane (void for none)
%
% OUTPUT :
%       - out_structure : the final structure with all the data
%                   The data structure is also saved in
%                   'fnc\data_reader\extracted_data\'output_str'.mat'
%warning('off');
disp('==============================')
disp('-- csv_DATEX_reader_v4 ')


try
    path=strcat(pwd,'\fnc\extracted_data\');
    min_freq = 6; % Set min_freq

    %% Load data
    % Data obtained with the "volledig" (full) structure
    input_str=strcat(pwd,'\fnc\data_reader\traffic_data\',input_str);
    filename = [input_str,'.csv'];
    fprintf('1) Using data in: %s \n',filename)

    % import in a cell
    import_raw = importdata(filename, ';');
    cell_raw = import_raw.textdata;
    index=1;
    for i = 1:size(cell_raw,2)
        if(isempty(cell_raw{2,i}))
            cell_raw(2:end,i) = num2cell(import_raw.data(:,index));
            index=index+1;
        end
    end

    % create an empty structure
    data = struct();

    disp('2) Reading data... ')
    %% Reading data: Header
    %The first row is the header, thus it is used to create the structure
    %fields
    header = cell_raw(1,:);

    % create all the fields
    for j = 1 : length(header)
        field_name = char(header(1,j));
        data.(field_name) = [];
    end

    %% Reading data: Data rows
    % fill all the fields row by row
    data_field_name = fieldnames(data);
    for i = 1:numel(data_field_name)
        % get the data field
        field_i_name = char(data_field_name(i));
        % select the data in the row
        data.(field_i_name) = [cell_raw(2:end, i)];
    end

    %% Find sensor names
    sensors_raw = unique(data.naam_meetlocatie_mst);
    for i = 1:length(sensors_raw)
        sensors_id(i) = erase(extractAfter(string(sensors_raw(i)), 8), 'ra');
    end

    %% Extract useful data
    % extract from the whole data only the values that interest us
    % find the indeces associated with the different sensors

    disp('3) Extracting useful data... ')
    sensor(length(sensors_id)) = struct(); %preallocate space for speed-up
    parfor j = 1:length(sensors_id)
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
        sensor(j).time_sample = str2double(data.gebruikte_minuten_intensiteit(sensor_index));
        sensor(j).ending_s_time = data.eind_meetperiode(sensor_index);
        sensor(j).starting_s_time = data.start_meetperiode(sensor_index);
        sensor(j).latitude = data.start_locatie_latitude(sensor_index);
        sensor(j).longitude = data.start_locatie_longitude(sensor_index);
        sensor(j).n_lanes = str2double(data.totaal_aantal_rijstroken(sensor_index));
    end

    %% Divide data between main lanes and service station lanes
    disp('4) Splitting main lanes and service station lanes... ')
    lane_ss = opt.laneSS;
    kk=1;
    jj=1;
    flag=false;
    flag2=false;
    flag3=false;
    for j = 1:length(sensors_id)
        k=1;
        p=1;
        for i = 1 : length(sensor(j).lane)
            app = sensor(j).lane(i);
            if(strcmp(lane_ss, app))
                flag3=true;
                sensor_ss(kk).id = sensor(j).id;
                sensor_ss(kk).vehicle_number(k) = sensor(j).veh_number(i);
                sensor_ss(kk).vehicle_speed(k) = sensor(j).veh_avg_speed(i);
                if(sensor_ss(kk).vehicle_speed(k)<0)
                    sensor_ss(kk).vehicle_speed(k)=0;
                end
                sensor_ss(kk).time_sample(k) = sensor(j).time_sample(i);
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
                sensor_main(jj).time_sample(p) = sensor(j).time_sample(i);
                sensor_main(jj).ending_s_time(p) = sensor(j).ending_s_time(i);
                sensor_main(jj).starting_s_time(p) = sensor(j).starting_s_time(i);
                sensor_main(jj).latitude(p) =  sensor(j).latitude(i);
                sensor_main(jj).longitude(p) = sensor(j).longitude(i);
                sensor_main(jj).lane(p) =  sensor(j).lane(i);
                sensor_main(jj).n_lanes(p) = sensor(j).n_lanes(i);
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
    clear sensor data
    if flag3
        save_file = [path, 'sensor_ss_no_interp','.mat'];
        save(save_file,'sensor_ss')
    end
    %% Main lanes merged in one single measure per time interval
    disp('5) Joining main lanes... ')
    sensor_sum(length(sensors_id)) = struct(); %preallocate space for speed-up
    for i = 1:length(sensors_id)
        k = 1;

        ff= false;
        if flag3
            for pp = 1:length(sensor_ss)
                if (strcmp(sensors_id(i), sensor_ss(pp).id))
                    ff = true;
                    break
                end
            end
        end
        if (ff)
            n_lanes = sensor_main(i).n_lanes(1)-1;
        else
            n_lanes = sensor_main(i).n_lanes(1);
        end

        for j = 1:n_lanes:length(sensor_main(i).starting_s_time)
            total_veh=0;
            for z=j:(j+n_lanes-1)
                total_veh = total_veh + sensor_main(i).veh_number(z);
            end
            if total_veh == 0
                w_avg_speed = 0;
            else
                total_speed = 0;
                for z=j:(j+n_lanes-1)
                    total_speed = total_speed + sensor_main(i).veh_number(z)*sensor_main(i).veh_avg_speed(z);
                end
                w_avg_speed = total_speed/total_veh;
            end
            total_veh = total_veh/sensor_main(i).n_lanes(j);
            sensor_sum(i).id = sensor_main(i).id;
            sensor_sum(i).latitude(k) = sensor_main(i).latitude(j);
            sensor_sum(i).longitude(k) = sensor_main(i).longitude(j);
            sensor_sum(i).vehicle_number(k) = total_veh;
            sensor_sum(i).vehicle_speed(k) = w_avg_speed;
            sensor_sum(i).ending_time(k) = sensor_main(i).ending_s_time(j);
            sensor_sum(i).starting_time(k) = sensor_main(i).starting_s_time(j);
            sensor_sum(i).sample_time(k) = sensor_main(i).time_sample(j);
            k = k+1;
        end
    end
    save_file = [path, 'sensor_sum_no_interp','.mat'];
    save(save_file,'sensor_sum')

    %% Interpolate the data
    disp('6) Reshaping and interpolating the data... ')
    % if the minimum frequency is higher than the data sample rate, we
    % interpolate the data to obtain the desired frequency

    if ~isempty(min_freq) && sensor_sum(1).sample_time(1) > 1/min_freq
        parfor k = 1: length(sensors_id)
            xx = 1:length(sensor_sum(k).vehicle_number);
            yy = 1:1/min_freq:length(sensor_sum(k).vehicle_number);

            veh_interp = interp1(xx,sensor_sum(k).vehicle_number,yy); %interpolate the veh flow
            speed_interp = interp1(xx,sensor_sum(k).vehicle_speed,yy); %interpolate the veh speed

            % sample time in [h]: from the website we get the data in
            % minutes, hence we need to divide by 60 to achieve the desired
            % measure
            sensor_sum(k).sample_time = sensor_sum(k).sample_time(1)/min_freq*ones(1,length(yy))*(1/60);

            for i=1:min_freq-1
                %extend the interpolated array to match the correct dimensions
                veh_interp = [veh_interp veh_interp(end)];
                speed_interp = [speed_interp speed_interp(end)];
                sensor_sum(k).sample_time = [sensor_sum(k).sample_time sensor_sum(k).sample_time(end)];
            end

            veh_interp = round(veh_interp);

            % assign the new values
            sensor_sum(k).vehicle_speed = speed_interp;
            sensor_sum(k).vehicle_number = veh_interp.*sensor_sum(k).sample_time;

            % extend the other fields in "sensor_sum"
            sensor_sum(k).starting_time = repelem(sensor_sum(k).starting_time, 1,min_freq);
            sensor_sum(k).ending_time = repelem(sensor_sum(k).ending_time, 1,min_freq);
            sensor_sum(k).latitude = repelem(sensor_sum(k).latitude, 1, min_freq);
            sensor_sum(k).longitude = repelem(sensor_sum(k).longitude, 1, min_freq);
        end
    end

    %% Compute the  fundamental diagram
    parfor j = 1:length(sensors_id)
        flow = sensor_sum(j).vehicle_number./sensor_sum(j).sample_time; % [veh/h]
        density = flow./sensor_sum(j).vehicle_speed;
        sensor_sum(j).flow = flow;
        sensor_sum(j).density = density;
    end
    %     f1 = figure;
    %     scatter(sensor_sum(10).density,sensor_sum(10).vehicle_speed,[],'filled')
    %     hold on
    %     grid on
    %     f1.WindowState = 'maximized';
    %     ax = gca();
    %     font_sz = 25;
    %     ax.XAxis.FontSize = font_sz; ax.XAxis.TickLabelInterpreter = 'latex';
    %     ax.YAxis.FontSize = font_sz; ax.YAxis.TickLabelInterpreter = 'latex';
    %     ax.XAxis.Label.String = '$\rho[veh/km]$'; ax.XAxis.Label.FontSize = font_sz;
    %     ax.XAxis.Label.Interpreter = 'latex';
    %     ax.YAxis.Label.String = '$speed[km/h]$'; ax.YAxis.Label.FontSize = font_sz;
    %     ax.YAxis.Label.Interpreter = 'latex';
    %     exportgraphics(f1,['density_speed.pdf'],...
    %                    'BackgroundColor','none');
    %     exportgraphics(f1,['density_speed.eps'],...
    %                    'BackgroundColor','none');

    %% Plot
    if(opt.display>0)
        last_fig_num = get(gcf,'Number');
        n_row = 3; N = size(sensor_sum,2);
        for n = 1 : N
            % % % % % % % %
            figure(last_fig_num+1)
            subplot(n_row,ceil((N)/n_row),n)
            bar(sensor_sum(n).vehicle_number)
            ax = gca;
            ax.YLim = [0,8];
            title_str1 = ['# vehicles (10s)-sens',sensor_sum(n).id];
            title(title_str1)
            grid on

            % % % % % % % %
            figure(last_fig_num+2)
            subplot(n_row,ceil((N)/n_row),n)
            bar(sensor_sum(n).vehicle_speed)
            grid on
            ax = gca;
            title_str2 = ['Avg. speed-sens',sensor_sum(n).id];
            title(title_str2)
            ax.YLim = [0,150];

            % % % % % % % %
            figure(last_fig_num+3)
            subplot(n_row,ceil((N)/n_row),n)
            scatter(sensor_sum(n).density,sensor_sum(n).flow)
            ax = gca;
            ax.YLim = [0,2500];
            grid on
            title_str3 = ['Fundamental diagram-sens',sensor_sum(n).id];
            title(title_str3)
        end
    end
    %% Save the file and assign the output
    save_file = [path, 'data','.mat'];
    save(save_file,'sensor_sum')

catch ME
    keyboard
    rethrow(ME)
end
end