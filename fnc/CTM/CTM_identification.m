function [CTM_param,phi_1] = CTM_identification(data,opt)
%CTM_identification : identify the parameters of the CTM model from the raw
%data extracted 
% INPUT:    - data : the structure with all the data about the traffic and
%                    cells               
%           - opt.disp {0,1} if we want to plot or not the figures
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
disp('==============================')

try
% number of the last figure
last_fig_num = get(gcf,'Number');    
% the number of sensors 
N = numel(data);
% sensors are placed at the interfaces of the cells, so
N_cell = N-1;
CTM_param.N = N_cell;

% speed thrashold
if length(opt.speed_th) == N_cell
    speed_th = opt.speed_th;
else
    error('ERROR : wrong dimension of "data.speed_th" ') 
end

if length(opt.coeff_quantile) ~= N_cell
    error('ERROR : wrong dimension of "data.coeff_quantile" ') 
end

%% Initialize the CT_param structure
CTM_param.v_bar = zeros(N_cell,1);
CTM_param.w = zeros(N_cell,1);
CTM_param.q_max = zeros(N_cell,1);
CTM_param.rho_max = zeros(N_cell,1);
CTM_param.rho_real = zeros(length(data(1).flow),N_cell);
CTM_param.len = zeros(N_cell,1);

for n = 1:N_cell
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
    position_s1 = str2num(cell2mat(data(n).position(1)));
    % position sensor n
    position_s2 = str2num(cell2mat(data(n+1).position(1)));
    CTM_param.len(n) = 10^(-3)*(abs(position_s2-position_s1)); % [km]
    CTM_param.T(n) = data(1).sample_time(1); % [h]
    
    %% Flow
    % In/out flows from the cell
    flow_IN = data(n).flow;
    if n == 1
        phi_1 = flow_IN;
    end
    flow_OUT = data(n+1).flow;
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
    ind_red_del = []; tol = 50000;
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
            fff = figure(last_fig_num+3);
            subplot(n_row,ceil((N-1)/n_row),n)
            yy = 1: length(density);
            scatter(yy,density,[],color,'filled')
            grid on
        end
end
%% Check if the frequency 
% the frequency has to be high enough
% check if it is needed 
if sum(CTM_param.T'.*CTM_param.v_bar./CTM_param.len>1)>0
    disp('==============================')
    disp('WARN: the frequency is too low!')
    disp(' - this might lead to a negative desity -')
    disp([' - to solve choose: T < ',num2str(min(CTM_param.len./CTM_param.v_bar)),' [h]   -']);
    disp('==============================')
end

catch ME
    keyboard
    rethrow(ME)
end
end
