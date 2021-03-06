function  out_structure = csv_DATEX_reader(input_str,output_str,opt,extra)
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

    %% Check number of onputs
    if nargin() < 4
        min_freq = [];
        sensors_id =  {'132','137','142','147'};
    else
        min_freq = extra.min_freq;
        sensors_id = extra.sensor_id;
        % check if the frequency id an integer or not
        if ~mod(min_freq,1) == 0
           error('ERROR : "min_freq" has to be an integer') 
        end
    end
    
    %% Check if the user wants to load other data
    % load a previously extracted data structure and immediately send it as
    % an output of the function
    if opt.load 
        file_to_load = ['fnc\data_reader\extracted_data\',...
                        output_str,'.mat'];
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
        cell_raw = importdata(filename);
        % create an empty structure
        data = struct();

        %% Header
        % the first row is the header, 
        % thus it is used to create the structure field 
        header = char(cell_raw(1));
        % find the commas in the field
        hd_comma = find(header == ','); 
        % num of commas 
        num_hd_commas = length(hd_comma);
        % if there is no comma at the end of the header then last elemnt in
        % hd_comma, we add as final index the size of the header
        hd_comma_end = [];
        if hd_comma ~= length(header)
           hd_comma_end = length(header)+1; 
        end
        hd_comma = [0, hd_comma, hd_comma_end];

        % create all the fields
        for j = 1: length(hd_comma)-1
            field_name = header(hd_comma(j)+1:hd_comma(j+1)-1);
            if ~isnan(str2double(field_name(1)))
                % to avoid error due to first el a number
                field_name = ['a', field_name];  
            end
            data.(field_name) = []; 
            % setfield(data,field_name,[]);
        end

        %% Fill the data fields
        if opt.verbatim
            fprintf('2) Construsction \n')
        end
        % complete all the fields with every row

        for k = 2:length(cell_raw)
            if mod(k,5000) == 0 && opt.verbatim
                fprintf('   -Row %d/%d \n',k,length(cell_raw))
            end
            % find the commas in the cell raw
            row = char(cell_raw(k));
            row_comma = find(row == ',');
            % add the 0 and the final element if necessary
            row_comma_end = [];
            if row_comma ~= length(row_comma)
               row_comma_end = length(row_comma)+1;
            end
            row_comma = [0, row_comma, row_comma_end];
            % 
            data_field_name = fieldnames(data);
            for i = 1:numel(data_field_name)
                % get the data field
                field_i_name = char(data_field_name(i)); 
                tmp = data.(field_i_name);
                % select the data in the row
                tmp2 = {row(row_comma(i)+1:row_comma(i+1)-1)};
                data.(field_i_name) = [tmp,tmp2];
            end
        end
        if opt.verbatim
            fprintf('3) Created struct : %s \n','data') 
        end
        %% Extract useful data
        % extract from the whole data only the ones that interest us
        % find different inde associated to the different sensors

        % sensors_id = {'132','137','142','147'};

        % find in the data the ones associated to a particular sensor
        for j = 1:length(sensors_id)
            sensor_name = char(sensors_id(j));
            check_sensor = strfind(data.naam_meetlocatie_mst,sensor_name);    
            tmp_sensor_pos = zeros(length(check_sensor),1);
            for i = 1:length(check_sensor)
                 tmp_sensor_pos(i) = ~isempty(cell2mat(check_sensor(i)));
            end
            sensor(j).data_index = tmp_sensor_pos;
        end
        %%
        for j = 1:length(sensors_id)
            sensor(j).id = char(sensors_id(j)); % transform to string
            ind = logical(sensor(j).data_index);
            sensor(j).veh_number = data.gem_intensiteit(ind); % the intensity of vehicles over one h
            sensor(j).veh_avg_speed = data.gem_snelheid(ind);
            %notice that the sample time is not always consistent
            % even though the most cases it is
            sensor(j).time_sample = data.gebruikte_minuten_intensiteit(ind);
            sensor(j).ending_s_time = data.eind_meetperiode(ind);
            sensor(j).starting_s_time = data.start_meetperiode(ind);
            sensor(j).location = data.afstand_tot_vild_locatie(ind);
            asd1 = zeros(sum(ind),1); asd2 = asd1; asd3 = asd1;
            for t = 1:sum(ind)
                asd1(t) = str2double(sensor(j).veh_number(t));
                asd2(t) = str2double(sensor(j).time_sample(t));        
                asd3(t) = str2double(sensor(j).veh_avg_speed(t));
            end
            sensor(j).veh_number = asd1;
            sensor(j).time_sample = asd2; 
            sensor(j).veh_avg_speed = asd3;    
        end
        if opt.verbatim
            fprintf('4) Organize in struct sensor : %s \n', 'sensor')
        end
        % collect all the data that are measured in the same tiem interval, they
        % are assumed consecutive
        if opt.verbatim
            fprintf('5) Reshape the data \n')
        end
        new_ind = 1;
        for k = 1:length(sensors_id)
            % check corrupted data, hence zeros that should not be there
            % and due to sensor failure
            for jj = 2 : length(sensor(k).veh_number)-1
                if sensor(k).veh_number(jj) == 0
                   % check if the previous one and one after are zero or not
                   if sensor(k).veh_number(jj-1) ~= 0 || sensor(k).veh_number(jj+1) ~= 0
                       sensor(k).veh_number(jj) = round((sensor(k).veh_number(jj-1)+sensor(k).veh_number(jj+1))/2);
                       sensor(k).veh_avg_speed(jj) = (sensor(k).veh_avg_speed(jj-1)+sensor(k).veh_avg_speed(jj+1))/2;
                   end
                end
                
            end
            
            % create some temporary structure  fields
            sensor(k).ending_time = sensor(k).ending_s_time;
            sensor(k).starting_time = sensor(k).ending_s_time;
            sensor(k).vehicle_number = sensor(k).veh_number;
            sensor(k).vehicle_speed = sensor(k).veh_avg_speed; 
            sensor(k).sample_time = sensor(k).time_sample;
            sensor(k).position = sensor(k).location;
            time_ind = '';
            for jj = 1 : length(sensor(k).starting_s_time)
                if strcmp(time_ind,sensor(k).ending_s_time(jj))
                    % it means that ti is the same time 
                    % create the weighted avg speed
                    total_veh = sensor(k).vehicle_number(new_ind)+sensor(k).veh_number(jj);
                    v1 = sensor(k).vehicle_speed(new_ind) ;
                    v2 = sensor(k).veh_avg_speed(jj);
                    % add the number of vehicles
                    sensor(k).vehicle_number(new_ind) = total_veh;
                    sensor(k).vehicle_speed(new_ind) = v1*sensor(k).vehicle_number(new_ind)/total_veh ...
                             +v2*sensor(k).veh_number(jj)/total_veh;
                else
                    % there sample period that starts
                    % save the starting and ending time
                    sensor(k).ending_time(new_ind) = sensor(k).ending_s_time(jj);
                    sensor(k).starting_time(new_ind) = sensor(k).starting_s_time(jj);
                    sensor(k).vehicle_number(new_ind) = sensor(k).veh_number(jj);
                    sensor(k).vehicle_speed(new_ind) = sensor(k).veh_avg_speed(jj);
                    sensor(k).sample_time(new_ind) = sensor(k).time_sample(jj); 
                    sensor(k).position(new_ind) = sensor(k).location(jj);
                    new_ind = new_ind+1;
                end
                time_ind = sensor(k).ending_s_time(jj);
            end
            sensor(k).ending_time(new_ind:end) = [];
            sensor(k).starting_time(new_ind:end) = [];
            sensor(k).vehicle_number(new_ind:end) = [];
            sensor(k).vehicle_speed(new_ind:end) = [];
            sensor(k).sample_time(new_ind:end) = []; 
            sensor(k).position(new_ind:end) = [];
            new_ind = 1;
            % make all column vectors
            sensor(k).ending_time = sensor(k).ending_time';
            sensor(k).starting_time = sensor(k).starting_time';
            sensor(k).position = sensor(k).position';
            
        end
        % remove the fields not used naymore
        sensor = rmfield(sensor,'ending_s_time');
        sensor = rmfield(sensor,'starting_s_time');
        sensor = rmfield(sensor,'time_sample');
        sensor = rmfield(sensor,'veh_number');
        sensor = rmfield(sensor,'veh_avg_speed');
        sensor = rmfield(sensor,'data_index');
        sensor = rmfield(sensor,'location');
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
                % [WRONG in v0.1:  We must divide by the frequency in order to obtain how
                % many veh travel through the new interval of time
                % interpolated_number_vv =
                % floor(interpolated_number_vv./min_freq); ]
                
                % This is  the correct one because "interpolated_number_vv"
                % is already in [veh/h]
                interpolated_number_vv = floor(interpolated_number_vv./min_freq);
                % interpolate the velocity
                speed_v = sensor(k).vehicle_speed;
                interpolated_speed_vv = interp1(xx,speed_v,yy);
                % assign the new values
                sensor(k).vehicle_number = interpolated_number_vv;
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
            end
        end

        %% Compute the  fundamental diagram
        % we already have the flow, we just need the density
        for j = 1:length(sensors_id)
            % [WRONG in v0.1: flow = vehicle_num*(sample_time/1h)
            % flow = sensor(j).vehicle_number./sensor(j).sample_time;%.*60;
            % % [veh/h]]
            
            % flow = vehicle_num -> it is already the flow cause it is the
            % number every hour already.
            flow = sensor(j).vehicle_number; % [veh/h]
            % density = flow*speed
            density = flow./sensor(j).vehicle_speed;
            sensor(j).flow = flow;
            sensor(j).density = density;
        end
        % assign the output
        out_structure = sensor;
        %% Save the file
        save_file = ['fnc\data_reader\extracted_data\',...
                        output_str,'.mat'];
        save(save_file,'sensor')
        if opt.verbatim
            fprintf('6) Save the data in %s\n',save_file)
            disp('==============================')
        end
    end
    %% Plot 
    try
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
    n_row = 2; N = size(sensor,2);
    for n = 1 : N
        % % % % % % % %
        figure(last_fig_num+1)
        subplot(n_row,ceil((N-1)/n_row),n)
        bar(sensor(n).vehicle_number)
        title_str1 = ['# vehicles (Sens. ',num2str(n),')'];
        title(title_str1)
        grid on
        
        % % % % % % % %
        figure(last_fig_num+2)
        subplot(n_row,ceil((N-1)/n_row),n)
        bar(sensor(n).vehicle_speed)
        grid on
        ax = gca;
        title_str2 = ['Avg. speed (Sens.',num2str(n),')'];
        title(title_str2)
        ax.YLim = [0,150];
        
        % % % % % % % %
        figure(last_fig_num+3)
        subplot(n_row,ceil((N-1)/n_row),n)
        scatter(sensor(n).density,sensor(n).flow)
        grid on
        title_str3 = ['Fundamental diagram (Sens. ',num2str(n),')'];
        title(title_str3)
    end
end






