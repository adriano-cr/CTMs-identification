function [] = CTMs_parameters(opt)
%CTMs_parameters : identify the parameters of the CTMs model from the raw
% data extracted
% INPUT:    - opt.disp {0,1} if we want to plot or not the figures
%           - opt.id_sensor_input
%           - opt.id_sensor_output

disp('==============================')
disp('-- CTMs parameters identification ')

try

    path=strcat(pwd,'\fnc\extracted_data\sensor_sum_no_interp.mat');
    aa = load(path, '*');
    sensor_sum = aa.sensor_sum;
    clear aa;
    
    %% Flows and beta estimation
    disp('1) Flows and beta estimation... ')
    id_sensor_input = opt.id_sensor_input;
    id_sensor_output = opt.id_sensor_output;

    for i = 1:length(sensor_sum)
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

    flow_in_clean = filloutliers(flow_in,"linear","percentiles",[25 99]);
    flow_out_clean = filloutliers(flow_out,"pchip","percentiles",[10 85]);

    x_flow = (linspace(0,24,length(flow_out_clean)))';
    [f_in_poly,gof] = fit(x_flow,flow_in_clean','poly6','Robust', 'Bisquare');
    [f_in_fou,gof] = fit(x_flow,flow_in_clean','fourier5','Robust', 'Bisquare');
    [f_out_poly,gof] = fit(x_flow,flow_out_clean','poly6','Robust', 'Bisquare');
    [f_out_fou,gof] = fit(x_flow,flow_out_clean','fourier5','Robust', 'Bisquare');

    beta = flow_in./(sensor_sum(index_sensor_input).vehicle_number);
    for i=1:length(beta)
        if(isnan(beta(i)))
            beta(i)=-1;
        end
    end

    x_beta=linspace(0,24,length(beta));
    beta_outliers=excludedata(x_beta,beta,'range',[0 1]);
    [f_beta,gof] = fit(x_beta',beta','poly2','Robust', 'Bisquare','Exclude', beta_outliers);



    %% Station occupancy estimation
    disp('2) Occupancy estimation... ')
    y_f_out_fou=f_out_fou(x_flow);
    y_f_in_fou=f_in_fou(x_flow);
    y_f_out_poly=f_out_poly(x_flow);
    y_f_in_poly=f_in_poly(x_flow);

    io_poly =  y_f_in_poly-y_f_out_poly;
    io_fou = y_f_in_fou-y_f_out_fou;
    occupancy_poly=[];
    occupancy_fou=[];

    window=60; % in minutes
    for i=1:window:length(y_f_in_poly)-window+1
        input=0;
        output=0;
        for j=i:i+window-1
            input = input + y_f_in_poly(j);
            output = output + abs(y_f_out_poly(j));
        end
        occupancy_poly = [occupancy_poly input-output];
        
        input=0;
        output=0;
        for j=i:i+window-1
            input = input + y_f_in_fou(j);
            output = output + abs(y_f_out_fou(j));
        end
        occupancy_fou = [occupancy_fou input-output];
    end

    %% Station delta estimation
    disp('3) Delta estimation... ')

    piecewise_out= ppcreate(x_flow,y_f_out_fou, 'pchip');
    piecewise_out= ppcreate(piecewise_out,'cut', [5.75 24]);
    inv_piecewise_out = ppcreate(piecewise_out,'inv');

    piecewise_in= ppcreate(x_flow,y_f_in_fou, 'pchip');
    piecewise_in= ppcreate(piecewise_in,'cut', [5.75 24]);
    inv_piecewise_in = ppcreate(piecewise_in,'inv');

    delay_poly=zeros(10,1);
    delay_fou=zeros(10,1);
    plot_tmp=[];
    plot_tmp2=[];
    ijk=1;

    for k=30:85
        %Fourier
        %inizialize variables
        x_in=[];
        x_out=[];
        time_input=[];
        time_output=[];
        [x_out,~]=inv_piecewise_out(k);
        [x_in,~]=inv_piecewise_in(k);

        if(length(x_out)>length(x_in))
            for i=1:length(x_in)
                for j=1:length(x_out)
                    if(x_out(j)>x_in(i))
                        time_output = [time_output x_out(j)];
                        time_input = [time_input x_in(i)];
                        break
                    end
                end
            end
        else 
            for i=length(x_out):-1:1
                for j=length(x_in):-1:1
                    if(x_out(i)>x_in(j))
                        time_input = [time_input x_in(j)];
                        time_output = [time_output x_out(i)];
                        break
                    end
                end
            end

        end
        app = 60.*(time_output-time_input);
        %copy delay to final matrix
        for line=1:length(app)
            delay_fou(line,ijk) = app(line);
            plot_tmp = [plot_tmp; app(line) (time_output(line)+time_input(line))/2];
        end

        %Polynomial
        %inizialize variables
        x_in=[];
        x_out=[];
        time_input=[];
        time_output=[];
        x_out_clean=[];
        x_in_clean=[];
        
        
        coefficientValues_f_in_poly = coeffvalues(f_in_poly);
        syms f_in_poly_sym(p1,p2,p3,p4,p5,p6,p7,x)
        eq=p1*x^6+p2*x^5+p3*x^4+p4*x^3+p5*x^2+p6*x+p7;
        f_in_poly_sym(p1,p2,p3,p4,p5,p6,p7,x)=eq;
        f_in_poly_sym(x) = subs(f_in_poly_sym(p1,p2,p3,p4,p5,p6,p7,x), {p1, p2, p3, p4, p5, p6, p7}, coefficientValues_f_in_poly);

        coefficientValues_f_out_poly = coeffvalues(f_out_poly);
        syms f_out_poly_sym(p1,p2,p3,p4,p5,p6,p7,x)
        f_out_poly_sym(p1,p2,p3,p4,p5,p6,p7,x)=eq;
        f_out_poly_sym(x) = subs(f_out_poly_sym(p1,p2,p3,p4,p5,p6,p7,x), {p1, p2, p3, p4, p5, p6, p7}, coefficientValues_f_out_poly);
        % solve equations to find roots
        x_out = double(solve(f_out_poly_sym(x) == k,x));
        x_in = double(solve(f_in_poly_sym(x) == k,x));

        % selecting only usuful roots
        for i=1:length(x_out)
            if((imag(x_out(i))==0)&&(real(x_out(i))>=0)&&(real(x_out(i))<=24))
                x_out_clean = [x_out_clean real(x_out(i))];
            end
            if((imag(x_in(i))==0)&&(real(x_in(i))>=0)&&(real(x_in(i))<=24))
                x_in_clean = [x_in_clean real(x_in(i))];
            end
        end

        if(~isempty(x_out_clean))&&(~isempty(x_in_clean))
            %usuful roots exist 
            if(length(x_out_clean)>length(x_in_clean))   
                for i=1:length(x_in_clean)
                    for j=1:length(x_out_clean)
                        if(x_out_clean(j)>x_in_clean(i))
                            time_output = [time_output x_out_clean(j)];
                            time_input = [time_input x_in_clean(i)];
                            break
                        end
                    end
                end
            else
                for i=length(x_out_clean):-1:1
                    for j=length(x_in_clean):-1:1
                        if(x_out_clean(i)>x_in_clean(j))
                            time_input = [time_input x_in_clean(j)];
                            time_output = [time_output x_out_clean(i)];
                            break
                        end
                    end
                end
            end
            app = 60.*(time_output-time_input);
        else
            %no usuful solutions
            app=[];
        end
        %copy delay to final matrix
        for line=1:length(app)
            delay_poly(line,ijk) = app(line);
            plot_tmp2 = [plot_tmp2; app(line) (time_output(line)+time_input(line))/2];
        end
        ijk=ijk+1; % update counter for the final matrix
        fprintf("\t - y = %d done\n", k)
    end

    % computations for plot - fourier
    delay_fou_plot = [];
    for i=1:length(plot_tmp)
        if(plot_tmp(i,1)*plot_tmp(i,2)>0)
            row = [plot_tmp(i,2) plot_tmp(i,1)];
            delay_fou_plot = [delay_fou_plot; row];
        end
    end
    [delay_fou_plot(:,1), I] = sort(delay_fou_plot(:,1));
    delay_fou_plot(:,2) = delay_fou_plot(I,2);

    % computations for plot - poly
    delay_poly_plot = [];
    for i=1:length(plot_tmp2)
        if(plot_tmp2(i,1)*plot_tmp2(i,2)>0)
            row = [plot_tmp2(i,2) plot_tmp2(i,1)];
            delay_poly_plot = [delay_poly_plot; row];
        end
    end
    [delay_poly_plot(:,1), I] = sort(delay_poly_plot(:,1));
    delay_poly_plot(:,2) = delay_poly_plot(I,2);

    %% Plots
    if(opt.display>0)
        last_fig_num = get(gcf,'Number');
        figure(last_fig_num+1)
        scatter(x_flow',flow_in_clean,'x')
        hold on
        plot(f_in_poly,'red')
        plot(f_in_fou, 'black')
        %ylim([-50 300])
        grid on
        xlabel("hour")
        ylabel("flow in")
        legend('dati', 'polinomio', 'fourier')
        title('flow in ripulito');

        figure(last_fig_num+2)
        scatter(x_flow',flow_out_clean,'x')
        hold on
        plot(f_out_poly,'red')
        plot(f_out_fou, 'black')
        %ylim([-10 300])
        legend('dati', 'polinomio', 'fourier')
        grid on
        xlabel("hour")
        ylabel("flow out")
        title('flow out ripulito');

        figure(last_fig_num+3)
        plot(f_in_poly,'red')
        grid on
        hold on
        plot(f_out_poly,'blue')

        %ylim([0 125])
        xlabel("hour")
        ylabel("flow")
        legend('flow in', 'flow out')
        title('poly flow in vs out puliti');

        figure(last_fig_num+4)
        plot(f_in_fou,'red')
        grid on
        hold on
        plot(f_out_fou,'blue')
        %ylim([0 125])
        legend('flow in', 'flow out')
        xlabel("hour")
        ylabel("flow")
        title('fourier flow in vs out puliti');

        figure(last_fig_num+5)
        plot(x_flow,io_poly)
        hold on
        plot(fit(x_flow,io_poly,'poly1'))
        legend('I/O poly', 'fitted curve')
        grid on
        xlabel("hour")
        ylabel("I/O service station")
        title('I/O Polynomial');

        figure(last_fig_num+6)
        hold on
        grid on
        plot(x_flow,io_fou)
        plot(fit(x_flow,io_fou,'poly1'))
        legend('I/O fourier', 'fitted curve')
        xlabel("hour")
        ylabel("I/O service station")
        title('I/O Fourier');

        figure(last_fig_num+7)
        hold on
        grid on
        w=24/length(occupancy_poly);
        aa=[];
        for i=0:w:24-w
            aa=[aa i];
        end
        plot(aa,occupancy_poly)
        plot(aa,occupancy_fou)
        legend('occupancy poly', 'occupancy fou')
        xticks(aa)
        xlabel("hour")
        ylabel("occupancy service station")
        title("Occupancy")

        figure(last_fig_num+8)
        subplot(2,1,1)
        scatter(delay_fou_plot(:,1), delay_fou_plot(:,2))
        hold on
        plot(fit(delay_fou_plot(:,1),delay_fou_plot(:,2),'poly1','Normalize','on','Robust','Bisquare'))
        %ylim([0 300])
        grid on
        ylabel("min")
        title('delta estimated with Fourier model');
        subplot(2,1,2)
        scatter(delay_poly_plot(:,1), delay_poly_plot(:,2))
        hold on
        plot(fit(delay_poly_plot(:,1),delay_poly_plot(:,2),'poly1','Normalize','on','Robust','Bisquare'))
        %ylim([0 300])
        ylabel("min")
        grid on
        title('delta estimated with poly model');
        
        figure(last_fig_num+9)
        scatter(x_beta,beta)
        hold on
        plot(f_beta, 'red')
        ylim([0 1])
        grid on
        xlabel("hour")
        title('beta');
    end
catch ME
    keyboard
    rethrow(ME)
end
end