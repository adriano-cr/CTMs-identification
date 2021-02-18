function [CTM_param] = CTM_identification(data,opt)
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
close all
% the number of sensors 
N = numel(data);
% sensors are placed at the interfaces of the cells, so
N_cell = N-1;
% speed thrashold
speed_th = 95;

%% Initialize the CT_param structure
CTM_param.v_bar = zeros(N,1);
CTM_param.w = zeros(N,1);
CTM_param.q_max = zeros(N,1);
CTM_param.rho_max = zeros(N,1);

for n = 1:N
    %% Categorize data point
    % extract the velocity 
    vel = data(n).vehicle_speed;
    % auxiliary structure
    color = [];
    congestion_bin = [];
    int_to_delete = [];

    for i = 1 : length(vel)
        if vel(i) < speed_th 
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

    %% Density and flow extraction
    % extract density
    density = data(n).density; density(int_to_delete) = [];
    % extract flow
    flow = data(n).flow; flow(int_to_delete) = [];
    % density/flow in free flow
    flow_green = flow(not_congestion_bin);
    density_green = density(not_congestion_bin);
    % density/flow in congestion flow
    density_red = density(congestion_bin);
    flow_red= flow(congestion_bin);

    % delete outliers
    flag = 0; i =1; ind_green_del =[];
    ind_red_del = []; tol = 10^-4;
    while flag == 0
        if i <= length(flow_green)
            if flow_green(i)<tol || density_green(i)<tol
                ind_green_del = [ind_green_del, i];
            end
        end
        % %
        if i <= length(flow_red)
            if flow_red(i)<tol || density_red(i)<tol
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
    f_green = [density_green\flow_green 0];
    p_green = polyval(f_green,density_green);
    % the velocity parameters
    CTM_param.v_bar(n) = f_green(1);

    % fit the congestion case
    [f_red,~]=quantreg(density_red,flow_red,.95,1);
    p_red = polyval(f_red,density_red);
    CTM_param.w(n) = -f_red(1); % minus cause f_red is negative
    CTM_param.rho_max(n) = -f_red(2)/f_red(1);

    % max flow of the cells, obtained as the intersection of congestion case
    % and free flow one
    CTM_param.q_max(n) = CTM_param.v_bar(n)*f_red(2)/(CTM_param.v_bar(n)+CTM_param.w(n));


    % cells length
    % CTM_param.len(n) = data.

    %% Plots
    n_row = 2;
        if opt.disp 
            % plot of cell i
            figure(1)
            subplot(n_row,ceil(N/n_row),n)
            xx = 1:length(data(n).vehicle_speed);
            scatter(xx,data(n).vehicle_speed,[],color,'filled')
            figure(2)
            subplot(n_row,ceil(N/n_row),n)
            scatter(density,flow,[],color,'filled')
            hold on
            yline(CTM_param.q_max(n),'b','LineWidth',3)
            plot(density_green,p_green,'b','LineWidth',3)
            plot(density_red,p_red,'b','LineWidth',3)
        end
end

end

