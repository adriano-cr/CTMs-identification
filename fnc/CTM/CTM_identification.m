function [CTM_param,phi_1,last_phi] = CTM_identification(opt)
%CTM_identification : identify the parameters of the CTM model from the raw
%data extracted
% INPUT:    - opt.disp {0,1} if we want to plot or not the figures
%             - opt.tolerance
%             - opt.speed_th
%             - opt.coeff_quantile
% OUTPUT:   - CTM_param.len
%           - CTM_param.v_bar
%           - CTM_param.w
%           - CTM_param.q_max
%           - CTM_param.rho_max
%           - CTM_param.len
%           - CTM_param.T
%           - CTM_param.N
%           - CTM_param.rho_real

disp('==============================')
disp('-- CTM identification ')


try
    disp('1) Initialize CTM param data... ')
    path=strcat(pwd,'\fnc\extracted_data\data.mat');
    aa = load(path, '*');
    data = aa.sensor_sum;
    clear aa;
    % number of the last figure
    last_fig_num = get(gcf,'Number');
    % the number of sensors
    N = numel(data);
    % sensors are placed at the interfaces of the cells, so
    N_cell = N-1;
    CTM_param.N = N_cell;
    % speed threshold
    speed_th = opt.speed_th;

    %% Initialize the CT_param structure
    CTM_param.v_bar = zeros(N_cell,1);
    CTM_param.w = zeros(N_cell,1);
    CTM_param.q_max = zeros(N_cell,1);
    CTM_param.rho_max = zeros(N_cell,1);
    CTM_param.rho_real = zeros(length(data(1).flow),N_cell);
    CTM_param.len = zeros(N_cell,1);
    disp('2) Cell computations: ')
    for n = 1:N_cell
        fprintf('\tAnalyzing cell %d/%d... \n', round(n),N_cell)
        %% Cell n velocity
        % extract the velocity as mean of the in and out velocity
        vel_s1 = data(n).vehicle_speed; % velocity sensor 1
        vel_s2 = data(n+1).vehicle_speed; % velocity sensor 2
        vel = (vel_s1+vel_s2)/2;
        vel = max(zeros(size(vel_s1)),vel);
        % auxiliary structure
        color = [];
        congestion_bin = [];

        % Colors: red = congested, green = free-flow
        for i = 1 : length(vel)
            if vel(i) < speed_th(n)
                % congestion case
                % red color
                color = [color; 1 0 0 ];
                % save as 1
                congestion_bin = [congestion_bin; 1];
            else
                % uncongested case
                % green color
                color = [color; 0 1 0 ];
                % save as 0
                congestion_bin = [congestion_bin; 0];
            end
        end
        % transform vector in boolean
        congestion_bin = boolean(congestion_bin);
        not_congestion_bin = boolean(1-congestion_bin);

        %% Cell length & time
        % position sensor n
        lat1=(cell2mat(data(n).latitude(1)));
        long1=(cell2mat(data(n).longitude(1)));
        lat2=(cell2mat(data(n+1).latitude(1)));
        long2=(cell2mat(data(n+1).longitude(1)));
        [CTM_param.len(n),~,~] = haversine([lat1 long1], [lat2 long2]);
        CTM_param.T(n) = data(1).sample_time(1); % [h]

        %% Flow
        % In/out flows from the cell
        flow_IN = data(n).flow;
        if n == 1
            phi_1 = flow_IN;
        end
        flow_OUT = data(n+1).flow;
        if n == N_cell
            last_phi = flow_OUT;
        end

        % average flow
        flow_avg = (flow_IN+flow_OUT)/2;

        %% Density
        density = max(flow_avg./vel,zeros(size(flow_avg)) );
        % density/flow in free flow
        flow_green = flow_avg(not_congestion_bin);
        density_green = density(not_congestion_bin);
        % density/flow in congestion flow
        density_red = density(congestion_bin);
        flow_red= flow_avg(congestion_bin);

        % delete outliers
        flag = 0; i =1; ind_green_del =[];
        ind_red_del = []; tol = opt.tolerance;
        while flag == 0
            if i <= length(flow_green)
                if flow_green(i)<0 || density_green(i)<0
                    ind_green_del = [ind_green_del, i];
                end
            end
            % %
            if i <= length(flow_red)
                if flow_red(i)<tol %% || density_red(i)<tol
                    ind_red_del = [ind_red_del, i];
                end
            end

            if i > length(flow_red) && i > length(flow_green)
                flag = 1;
            end
            i = i +1;
        end
        flow_red(ind_red_del) = []; density_red(ind_red_del)=[];
        flow_green(ind_green_del) = []; density_green(ind_green_del)=[];

        %% Identify the CTM parameters
        % fit the free flow case
        f_green = [mean(density_green.\flow_green) 0];
        p_green = polyval(f_green,density_green);
        % the velocity parameters
        CTM_param.v_bar(n) = f_green(1);

        % fit the congestion case
        coeff_quantreg = opt.coeff_quantile;
        [f_red,~]= quantreg(density_red',flow_red',coeff_quantreg(n),1);

        p_red = polyval(f_red,density_red);

        CTM_param.w(n) = -f_red(1); % minus cause f_red is negative
        CTM_param.rho_max(n) = -f_red(2)/f_red(1);

        % max flow of the cells, obtained as the intersection of congestion case
        % and free flow one
        CTM_param.q_max(n) = CTM_param.v_bar(n)*f_red(2)/(CTM_param.v_bar(n)+CTM_param.w(n));

        % cells density
        CTM_param.rho_real(:,n) = density;

        %% Compute the supply of cell N+1
        % create a variable supply of the last cell that will be equal to
        % q_max^{N} if there is no congestion from the data and lower
        % otherwise
        if n == N_cell % hence the last cell
            supply_N_plus_computation;
            CTM_param.supply_N_plus = supply_N_plus;
        end
        %% Plots
        n_row = 2;
        if opt.disp
            % plot of cell i
            f = figure(last_fig_num+1);
            subplot(n_row,ceil((N-1)/n_row),n)
            xx = 1:length(vel);
            scatter(xx,vel,[],color,'filled')
            title_str3 = ['Velocity (Cell ',num2str(n),')'];
            title(title_str3)
            grid on

            ff = figure(last_fig_num+2);
            subplot(n_row,ceil((N-1)/n_row),n)
            scatter(density_green',flow_green',[],'g','filled')
            hold on
            scatter(density_red,flow_red,[],'r','filled')
            % scatter(density,flow_avg,[],color,'filled')
            hold on
            grid on
            plot([0,max(density)],...
                [CTM_param.q_max(n),CTM_param.q_max(n)],'b','LineWidth',3)
            plot(density_green,p_green,'y','LineWidth',3)
            plot(density_red,p_red,'b','LineWidth',3)
            title_str3 = ['(Cell ',num2str(n),')'];
            title(title_str3)

            fff = figure(last_fig_num+3);
            subplot(n_row,ceil((N-1)/n_row),n)
            yy = 1: length(density);
            scatter(yy,density,[],color,'filled')
            title_str3 = ['Density (Cell ',num2str(n),')'];
            title(title_str3)
            grid on
        end
    end

    % %% Nice plot
    % % Fundamental diagram sensor 1
    % %
    % f1 = figure;
    % scatter(density_green',flow_green',[],'g','filled')
    % hold on
    % grid on
    % scatter(density_red,flow_red,[],'r','filled')
    % plot([0,max(density)],...
    %         [CTM_param.q_max(1),CTM_param.q_max(1)],'c','LineWidth',3)
    % plot(density_green,p_green,'y','LineWidth',3)
    % plot(density_red,p_red,'b','LineWidth',3)
    % f1.WindowState = 'maximized';
    % ax = gca();
    % font_sz = 25;
    % ax.XAxis.FontSize = font_sz; ax.XAxis.TickLabelInterpreter = 'latex';
    % ax.YAxis.FontSize = font_sz; ax.YAxis.TickLabelInterpreter = 'latex';
    % ax.XAxis.Label.String = '$\rho_1$'; ax.XAxis.Label.FontSize = font_sz;
    % ax.XAxis.Label.Interpreter = 'latex';
    % ax.YAxis.Label.String = '$\phi_1$'; ax.YAxis.Label.FontSize = font_sz;
    % ax.YAxis.Label.Interpreter = 'latex';
    % %exportgraphics(f1,['figure\fundamental_diagram_1.pdf'],...
    %  %               'BackgroundColor','none');
    % %exportgraphics(f1,['figure\fundamental_diagram_1.eps'],...
    % %                'BackgroundColor','none');
    %
    % %
    % f2 = figure;
    % xx = 1:length(vel);scatter(xx,vel,[],color,'filled')
    % grid on
    % f2.WindowState = 'maximized';
    % ax = gca();
    % K = length(vel);
    % font_sz = 25;
    %
    % ax.XTick = [1 round(1/CTM_param.T(1)*[1:23]) K];
    % asdasd = ax.XTick+1; asdasd(1) = asdasd(1)-1; asdasd(end) = asdasd(end)-1;
    % hr_lables = CS_status.time_hr(asdasd);
    % hr_string = datestr(hr_lables,'HH:MM');
    % select only between '06:00' and '22:00'
    % hr_string = hr_string(4:end-1,:);
    % asdasd = asdasd(4:end-1);
    % instants = [8:2:20,21];
    % hr_string = hr_string(instants,:);
    % asdasd = asdasd(instants);
    % ax.XTick = asdasd;
    % hr_cell = cellstr(hr_string); ax.XTickLabel = hr_cell;
    % ax.XAxis.FontSize = font_sz; ax.XAxis.TickLabelInterpreter = 'latex';
    % ax.YAxis.FontSize = font_sz; ax.YAxis.TickLabelInterpreter = 'latex';
    % ax.XLim = [asdasd(1), asdasd(end)];
    % ax.XAxis.FontSize = font_sz; ax.XAxis.TickLabelInterpreter = 'latex';
    % ax.YAxis.FontSize = font_sz; ax.YAxis.TickLabelInterpreter = 'latex';
    % ax.XAxis.Label.String = 'time $[\textup h]$'; ax.XAxis.Label.FontSize = font_sz;
    % ax.XAxis.Label.Interpreter = 'latex';
    % ax.YAxis.Label.String = '$v_1\,[\textup{km/h}]$'; ax.YAxis.Label.FontSize = font_sz;
    % ax.YAxis.Label.Interpreter = 'latex';
    % % exportgraphics(f2,['figure\velocity_1.pdf'],...
    % %                 'BackgroundColor','none');
    % % exportgraphics(f2,['figure\velocity_1.eps'],...
    % %                 'BackgroundColor','none');

catch ME
    keyboard
    rethrow(ME)
end
end