function  out_structure = csv_DATEX_reader_v4(input_str,output_str,opt)
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
%       - opt.display {0,1} : plot or not the foundamental graphs associated
%       to the traffic
%        - opt.laneSS
%        - opt.id_sensor_input
%        - opt.id_sensor_output
% OUTPUT :
%       - out_structure : the final structure with all the data
%                   The data structure is also saved in
%                   'fnc\data_reader\extracted_data\'output_str'.mat'

path = '\traffic_data';
addpath(genpath([pwd,path]))
min_freq = 6; % Set min_freq
try
    %% Load data
    % Data obtained with the "volle" (full) structure
    filename = [input_str,'.csv'];

    disp('==============================')
    fprintf('1) Using data in: %s \n',filename)
    disp('==============================')

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
    %field
    header = cell_raw(1,:);

    % create all the fields
    for j = 1 : length(header)
        field_name = char(header(1,j));
        data.(field_name) = [];
    end

    %% Reading data: Data rows
    % complete all the fields with every row
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
    % extract from the whole data only the ones that interest us
    % find different inde associated to the different sensors

    disp('3) Extracting useful data... ')
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
        sensor(j).time_sample = str2double(data.gebruikte_minuten_intensiteit(sensor_index));
        sensor(j).ending_s_time = data.eind_meetperiode(sensor_index);
        sensor(j).starting_s_time = data.start_meetperiode(sensor_index);
        sensor(j).latitude = data.start_locatie_latitude(sensor_index);
        sensor(j).longitude = data.start_locatie_longitude(sensor_index);
        sensor(j).n_lanes = str2double(data.totaal_aantal_rijstroken(sensor_index));
    end

    %% check errors due sensors failure
    % Done manually by selecting a day in which there are no errors or
    % sensors failure

    %% divide data for main lanes and service station lanes
    disp('4) Splitting main lanes and service station lanes... ')
    lane_ss = opt.laneSS;
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


    %% join the main lanes in one single measure per time interval
    disp('5) Joining main lanes... ')
    sensor_sum(length(sensors_id)) = struct(); %preallocate space for speed-up
    for i = 1:length(sensors_id)
        k = 1;
        for j = 1:5:length(sensor_main(i).starting_s_time)
            veh1 = sensor_main(i).veh_number(j);
            veh2 = sensor_main(i).veh_number(j+1);
            veh3 = sensor_main(i).veh_number(j+2);
            veh4 = sensor_main(i).veh_number(j+3);
            veh5 = sensor_main(i).veh_number(j+4);

            total_veh = veh1 + veh2 + veh3 + veh4 + veh5;
            %total_veh = veh1 + veh2 + veh3;

            vel1 = sensor_main(i).veh_avg_speed(j);
            vel2 = sensor_main(i).veh_avg_speed(j+1);
            vel3 = sensor_main(i).veh_avg_speed(j+2);
            vel4 = sensor_main(i).veh_avg_speed(j+3);
            vel5 = sensor_main(i).veh_avg_speed(j+4);
            if total_veh == 0
                w_avg_speed = 0;
            else
                w_avg_speed = (vel1*veh1 + vel2*veh2 + vel3*veh3 + veh4*vel4 + veh5*vel5)/total_veh;
                %w_avg_speed = (vel1 * veh1 + vel2 * veh2 + vel3 * veh3)/total_veh;
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
    %% calculate CTM-s parameters
    disp('6) Reshaping and interpolating the data... ')
    id_sensor_input = opt.id_sensor_input;
    id_sensor_output = opt.id_sensor_output;


    for i = 1:length(sensors_id)
        if(sensor_sum(i).id == id_sensor_output)
            index_sensor_output=i;
        end
        if(sensor_sum(i).id == id_sensor_input)
            index_sensor_input=i;
        end
    end


    flow_in=[];
    flow_out=[];
    for i=1: length(sensor_sum(1).vehicle_number)
        flow_out = [flow_out (sensor_sum(index_sensor_output).vehicle_number(i) - sensor_sum(index_sensor_output-1).vehicle_number(i))];
        flow_in = [flow_in (sensor_sum(index_sensor_input).vehicle_number(i) - sensor_sum(index_sensor_input+1).vehicle_number(i))];

    end

    flow_in_clean = filloutliers(flow_in,"nearest","percentiles",[40 98]);
    flow_out_clean = filloutliers(flow_out,"nearest","percentiles",[20 85]);

    %delay = flow_out_clean - flow_in_clean;
    beta = flow_in./(sensor_sum(index_sensor_input).vehicle_number);

    for i=1:length(beta)
        if(isnan(beta(i)))
            beta(i)=-1;
        end
    end
   
    if(opt.display>0)
        x=linspace(1,24,length(flow_out_clean));
        f_in_poly = fit(x',flow_in_clean','poly6');
        f_out_poly = fit(x',flow_out_clean','poly6');
        f_in_fou = fit(x',flow_in_clean','fourier5');
        f_out_fou = fit(x',flow_out_clean','fourier5');
        x_beta=linspace(1,24,length(beta));
        beta_outliers=excludedata(x_beta,beta,'range',[0 1]);
        f_beta = fit(x_beta',beta','fourier2', 'Exclude', beta_outliers);
         
        figure(1)
        scatter(x,flow_in_clean,'x')
        hold on
        plot(f_in_poly,'red')
        plot(f_in_fou, 'black')
        grid on
        legend('dati', 'polinomio', 'fourier')
        title('flow in ripulito');

        figure(2)
        scatter(x,flow_out_clean,'x')
        hold on
        plot(f_out_poly,'red')
        plot(f_out_fou, 'black')
        legend('dati', 'polinomio', 'fourier')
        grid on
        title('flow out ripulito');

        figure(3)
        plot(f_in_poly,'red')
        grid on
        hold on
        plot(f_out_poly,'blue')
        legend('flow in', 'flow out')
        title('poly flow in vs out puliti');

        figure(4)
        plot(f_in_fou,'red')
        grid on
        hold on
        plot(f_out_fou,'blue')
        legend('flow in', 'flow out')
        title('fourier flow in vs out puliti');

        figure(5)
        scatter(x_beta,beta)
        hold on
        plot(f_beta, 'red')
        grid on
        title('beta');
    end
    %% Interpolate the data
    % if the minimum frequency is higher than the one of the
    % data we interpolate the data to attain the desired one
    if ~isempty(min_freq) && sensor_sum(1).sample_time(1) > 1/min_freq  
        for k = 1: length(sensors_id)
            xx = 1:length(sensor_sum(k).vehicle_number); 
            yy = 1:1/min_freq:length(sensor_sum(k).vehicle_number);
          
            veh_interp = interp1(xx,sensor_sum(k).vehicle_number,yy); %interpolate the veh flow
            speed_interp = interp1(xx,sensor_sum(k).vehicle_speed,yy); %interpolate the veh speed

            % sample time in [h], from the site we have the data in
            % minutes hence we have to multiply 1/60 to achieve
            sensor_sum(k).sample_time = sensor_sum(k).sample_time(1)/min_freq*ones(1,length(yy))*(1/60);

            for i=1:min_freq-1 
                %extend the array interpolated to correctly match the dimensions
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
    for j = 1:length(sensors_id)
        flow = sensor_sum(j).vehicle_number./sensor_sum(j).sample_time; % [veh/h]
        density = flow./sensor_sum(j).vehicle_speed;
        sensor_sum(j).flow = flow;
        sensor_sum(j).density = density;
    end
    %% Plot
    if(opt.display>0)
        traffic_data_plot(sensor_sum)
    end
    %% Save the file and assign the output
    save_file = [opt.path, output_str,'.mat'];
    save(save_file,'sensor')
    fprintf('7) Saving the data in %s...\n',save_file)
    disp('==============================')
    out_structure = sensor_sum;

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
    bar(sensor_sum(n).vehicle_number)
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
    grid on
    title_str3 = ['Fundamental diagram-sens',sensor_sum(n).id];
    title(title_str3)
end
end