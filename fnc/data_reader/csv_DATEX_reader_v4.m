function  out_structure = csv_DATEX_reader_v4(input_str,output_str,opt)
%% csv_DATEX_reader :
% Function that reads the raw traffic data file and then creates a
% structure as output with all the data nicely organized and ready to be
% elaborated. It also displays some plots of the data if required
% (the plots in this case are not automatically saved).
%
% INPUT :
%       - input_str: the string that define the input file where the raw
%       data are extracted
%       - output_str: the name of the file in which we want to save the
%       strucute that we compute
%       - opt.display {0,1}: 1 to plot the foundamental graphs associated
%       with the traffic data, 0 to omit
%        - opt.laneSS: name of service station lane (void for none)
%        - opt.id_sensor_input: code of sensor before service station
%        - opt.id_sensor_output: code of sensor after service station
% OUTPUT :
%       - out_structure : the final structure with all the data
%                   The data structure is also saved in
%                   'fnc\data_reader\extracted_data\'output_str'.mat'

path = '\traffic_data';
addpath(genpath([pwd,path]))
min_freq = 6; % Set min_freq
try
    %% Load data
    % Data obtained with the "volledig" (full) structure
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
    % extract from the whole data only the values that interest us
    % find different indeces associated with the different sensors

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

    %% Check for errors due to sensors failure
    % Done manually by selecting a day with no errors or sensors failure

    %% Divide data between main lanes and service station lanes
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


    %% Main lanes merged in one single measure per time interval
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
    
    %% Estimation of CTM-s parameters
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

    %% Station delta estimation

    x_flow=(linspace(0,24,length(flow_out_clean)))';
    f_in_poly = fit(x_flow,flow_in_clean','poly6');
    f_out_poly = fit(x_flow,flow_out_clean','poly6');
    f_in_fou = fit(x_flow,flow_in_clean','fourier5');
    f_out_fou = fit(x_flow,flow_out_clean','fourier5');

    y_f_out_fou=round(f_out_fou(x_flow),2);
    y_f_in_fou=round(f_in_fou(x_flow),2);

    delay_poly=zeros(10,36);
    delay_fou=zeros(10,36);
    plot_tmp=[];
    ijk=1;

    for k=45:80
        delay_value = k;
        diff_in=mod(y_f_in_fou,delay_value);
        diff_out=mod(y_f_out_fou,delay_value);

        maxval_in=max(abs(diff_in-delay_value));
        tol=0.5;
        idx_in=[];
        for i=1:length(diff_in)
            if((diff_in(i)<=maxval_in+tol)&&(diff_in(i)>=maxval_in-tol))
                idx_in=[idx_in i];
            end
        end

        idx_in2=[];
        for i=1:length(idx_in)-1

            if(idx_in(i)+1~=idx_in(i+1))
                idx_in2=[idx_in2 idx_in(i)];
            end
            if(i==length(idx_in)-1)
                if(idx_in2(end)+1~=idx_in(end))
                    idx_in2=[idx_in2 idx_in(end)];
                end
            end
        end


        maxval_out=max(abs(diff_out-delay_value));
        tol=0.5;
        idx_out=[];
        for i=1:length(diff_in)
            if((diff_out(i)<=maxval_out+tol)&&(diff_out(i)>=maxval_out-tol))
                idx_out=[idx_out i];
            end
        end
        idx_out2=[];
        for i=1:length(idx_out)-1

            if(idx_out(i)+1~=idx_out(i+1))
                idx_out2=[idx_out2 idx_out(i)];
            end
            if(i==length(idx_out)-1)
                if(idx_out2(end)+1~=idx_out(end))
                    idx_out2=[idx_out2 idx_out(end)];
                end
            end
        end


        time_input=[];
        time_output=[];
        if(length(idx_out2)<length(idx_in2))
            for i=length(idx_out2):-1:1
                time_output = [time_output x_flow(idx_out2(i))];
                for j=length(idx_in2):-1:1
                    if(idx_in2(j)<idx_out2(i))
                        time_input = [time_input x_flow(idx_in2(j))];
                        break
                    end
                end
            end
        else
            for i=1:length(idx_in2)
                time_input = [time_input x_flow(idx_in2(i))];
                for j=1:length(idx_out2)
                    if(idx_out2(j)>idx_in2(i))
                        time_output = [time_output x_flow(idx_out2(j))];
                        break
                    end
                end
            end
        end

        app = 60.*(time_output-time_input);
        for line=1:length(app)
            delay_fou(line,ijk) = app(line);
        end
        coefficientValues_f_in_poly = coeffvalues(f_in_poly);
        syms f_in_poly_sym(p1,p2,p3,p4,p5,p6,p7,x)
        f_in_poly_sym(p1,p2,p3,p4,p5,p6,p7,x)=p1*x^6+p2*x^5+p3*x^4+p4*x^3+p5*x^2+p6*x+p7;
        f_in_poly_sym(x) = subs(f_in_poly_sym(p1,p2,p3,p4,p5,p6,p7,x), {p1, p2, p3, p4, p5, p6, p7}, coefficientValues_f_in_poly);

        coefficientValues_f_out_poly = coeffvalues(f_out_poly);
        syms f_out_poly_sym(p1,p2,p3,p4,p5,p6,p7,x)
        f_out_poly_sym(p1,p2,p3,p4,p5,p6,p7,x)=p1*x^6+p2*x^5+p3*x^4+p4*x^3+p5*x^2+p6*x+p7;
        f_out_poly_sym(x) = subs(f_out_poly_sym(p1,p2,p3,p4,p5,p6,p7,x), {p1, p2, p3, p4, p5, p6, p7}, coefficientValues_f_out_poly);

        eqn_out =  f_out_poly_sym(x) == delay_value;
        eqn_in =  f_in_poly_sym(x) == delay_value;

        sols_out = double(solve(eqn_out,x));
        sols_in = double(solve(eqn_in,x));

        sols_out_clean=[];
        sols_in_clean=[];
        
        for i=1:length(sols_out)
            if((imag(sols_out(i))==0)&&(real(sols_out(i))>=0)&&(real(sols_out(i))<=24))
                sols_out_clean = [sols_out_clean real(sols_out(i))];
            end
            if((imag(sols_in(i))==0)&&(real(sols_in(i))>=0)&&(real(sols_in(i))<=24))
                sols_in_clean = [sols_in_clean real(sols_in(i))];
            end
        end

        if(~isempty(sols_out_clean))&&(~isempty(sols_in_clean))
            app = 60.*(sols_out_clean-sols_in_clean);

        else
            app=[];
            fprintf("Useful solutions not available for y = %d\n", delay_value)
        end
        for line=1:length(app)
            delay_poly(line,ijk) = app(line);
            plot_tmp = [plot_tmp; app(line) (sols_out_clean+sols_in_clean)/2];
        end

        ijk=ijk+1;
    end
    
    delay_poly_plot = [];
    for i=1:length(plot_tmp)
        row1 = [plot_tmp(i,3) plot_tmp(i,1)];
        row2 = [plot_tmp(i,3) plot_tmp(i,2)];
        delay_poly_plot = [delay_poly_plot; row1; row2];
    end

    x_beta=linspace(0,24,length(beta));
    beta_outliers=excludedata(x_beta,beta,'range',[0 1]);
    f_beta = fit(x_beta',beta','fourier2', 'Exclude', beta_outliers);
    
    %% Plots
    if(opt.display>0)
        figure(1)
        scatter(x_flow',flow_in_clean,'x')
        hold on
        plot(f_in_poly,'red')
        plot(f_in_fou, 'black')
        xlim=[0 24];
        grid on
        legend('dati', 'polinomio', 'fourier')
        title('flow in ripulito');

        figure(2)
        scatter(x_flow',flow_out_clean,'x')
        hold on
        xlim=[0 24];
        plot(f_out_poly,'red')
        plot(f_out_fou, 'black')
        legend('dati', 'polinomio', 'fourier')
        grid on
        title('flow out ripulito');

        figure(3)
        plot(f_in_poly,'red')
        xlim=[0 24];
        grid on
        hold on
        plot(f_out_poly,'blue')
        legend('flow in', 'flow out')
        title('poly flow in vs out puliti');

        figure(4)
        plot(f_in_fou,'red')
        grid on
        hold on
        xlim=[0 24];
        plot(f_out_fou,'blue')
        legend('flow in', 'flow out')
        title('fourier flow in vs out puliti');

        figure(5)
        scatter(x_beta,beta)
        xlim=[0 24];
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