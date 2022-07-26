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
    last_fig_num = get(gcf,'Number');
    %% Estimation of CTM-s parameters
    disp('1) Flows estimation... ')
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

    flow_in_clean = filloutliers(flow_in,"nearest","percentiles",[40 98]);
    flow_out_clean = filloutliers(flow_out,"nearest","percentiles",[20 85]);
    beta = flow_in./(sensor_sum(index_sensor_input).vehicle_number);
    
    for i=1:length(beta)
        if(isnan(beta(i)))
            beta(i)=-1;
        end
    end

    %% Station delta estimation
    disp('2) Delta estimation... ')
    x_flow=(linspace(0,24,length(flow_out_clean)))';
    f_in_poly = fit(x_flow,flow_in_clean','poly6');
    f_out_poly = fit(x_flow,flow_out_clean','poly6');
    f_in_fou = fit(x_flow,flow_in_clean','fourier5');
    f_out_fou = fit(x_flow,flow_out_clean','fourier5');
    y_f_out_fou=f_out_fou(x_flow);
    y_f_in_fou=(f_in_fou(x_flow));

    piecewise_out= ppcreate(x_flow,y_f_out_fou, 'pchip');
    piecewise_out= ppcreate(piecewise_out,'cut', [5.75 24]);
    inv_piecewise_out = ppcreate(piecewise_out,'inv');

    piecewise_in= ppcreate(x_flow,y_f_in_fou, 'pchip');
    piecewise_in= ppcreate(piecewise_in,'cut', [5.75 24]);
    inv_piecewise_in = ppcreate(piecewise_in,'inv');

    delay_poly=zeros(10,36);
    delay_fou=zeros(10,36);
    plot_tmp=[];
    plot_tmp2=[];
    ijk=1;

    for k=20:80
        [x_out,~]=inv_piecewise_out(k);
        [x_in,~]=inv_piecewise_in(k);
        time_input=[];
        time_output=[];

        if(length(x_out)>length(x_in))
            %x_out longer
            for i=1:length(x_in)
                time_input = [time_input x_in(i)];
                for j=1:length(x_out)
                    if(x_out(j)>x_in(i))
                        time_output = [time_output x_out(j)];
                        break
                    end
                end
            end
        else
            %x_in longer
            for i=length(x_out):-1:1
                time_output = [time_output x_out(i)];
                for j=length(x_in):-1:1
                    if(x_out(i)>x_in(j))
                        time_input = [time_input x_in(j)];
                        break
                    end
                end
            end

        end
        app = 60.*(time_output-time_input);
        for line=1:length(app)
            delay_fou(line,ijk) = app(line);
            plot_tmp = [plot_tmp; app(line) (time_output(line)+time_input(line))/2];
        end

        coefficientValues_f_in_poly = coeffvalues(f_in_poly);
        syms f_in_poly_sym(p1,p2,p3,p4,p5,p6,p7,x)
        f_in_poly_sym(p1,p2,p3,p4,p5,p6,p7,x)=p1*x^6+p2*x^5+p3*x^4+p4*x^3+p5*x^2+p6*x+p7;
        f_in_poly_sym(x) = subs(f_in_poly_sym(p1,p2,p3,p4,p5,p6,p7,x), {p1, p2, p3, p4, p5, p6, p7}, coefficientValues_f_in_poly);

        coefficientValues_f_out_poly = coeffvalues(f_out_poly);
        syms f_out_poly_sym(p1,p2,p3,p4,p5,p6,p7,x)
        f_out_poly_sym(p1,p2,p3,p4,p5,p6,p7,x)=p1*x^6+p2*x^5+p3*x^4+p4*x^3+p5*x^2+p6*x+p7;
        f_out_poly_sym(x) = subs(f_out_poly_sym(p1,p2,p3,p4,p5,p6,p7,x), {p1, p2, p3, p4, p5, p6, p7}, coefficientValues_f_out_poly);

        eqn_out =  f_out_poly_sym(x) == k;
        eqn_in =  f_in_poly_sym(x) == k;

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
            fprintf("Useful solutions not available for y = %d\n", k)
        end
        for line=1:length(app)
            delay_poly(line,ijk) = app(line);
            plot_tmp2 = [plot_tmp2; app(line) (sols_out_clean(line)+sols_in_clean(line))/2];
        end

        ijk=ijk+1;
    end

    delay_fou_plot = [];
    for i=1:length(plot_tmp)
        if(plot_tmp(i,1)*plot_tmp(i,2)>0)
            row = [plot_tmp(i,2) plot_tmp(i,1)];
            delay_fou_plot = [delay_fou_plot; row];
        end
    end

    delay_poly_plot = [];
    for i=1:length(plot_tmp2)
        if(plot_tmp2(i,1)*plot_tmp2(i,2)>0)
            row = [plot_tmp2(i,2) plot_tmp2(i,1)];
            delay_poly_plot = [delay_poly_plot; row];
        end
    end

    [delay_fou_plot(:,1), I] = sort(delay_fou_plot(:,1));
    delay_fou_plot(:,2) = delay_fou_plot(I,2);

    [delay_poly_plot(:,1), I] = sort(delay_poly_plot(:,1));
    delay_poly_plot(:,2) = delay_poly_plot(I,2);

    x_beta=linspace(0,24,length(beta));
    beta_outliers=excludedata(x_beta,beta,'range',[0 1]);
    f_beta = fit(x_beta',beta','fourier2', 'Exclude', beta_outliers);

    %% Plots
    if(opt.display>0)
        figure(last_fig_num+1)
        scatter(x_flow',flow_in_clean,'x')
        hold on
        plot(f_in_poly,'red')
        plot(f_in_fou, 'black')
        ylim([-50 300])
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
        ylim([-50 300])
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
        ylim([0 100])
        xlabel("hour")
        ylabel("flow")
        legend('flow in', 'flow out')
        title('poly flow in vs out puliti');

        figure(last_fig_num+4)
        plot(f_in_fou,'red')
        grid on
        hold on
        plot(f_out_fou,'blue')
        ylim([0 100])
        legend('flow in', 'flow out')
        xlabel("hour")
        ylabel("flow")
        title('fourier flow in vs out puliti');

        figure(last_fig_num+5)
        scatter(x_beta,beta)
        hold on
        plot(f_beta, 'red')
        ylim([0 1])
        grid on
        xlabel("hour")
        title('beta');

        figure(last_fig_num+6)
        subplot(2,1,1)
        scatter(delay_fou_plot(:,1), delay_fou_plot(:,2))
        ylim([0 300])
        grid on
        ylabel("min")
        title('delta estimated with Fourier model');
        subplot(2,1,2)
        scatter(delay_poly_plot(:,1), delay_poly_plot(:,2))
        ylim([0 300])
        ylabel("min")
        grid on
        title('delta estimated with poly model');
    end
catch ME
    keyboard
    rethrow(ME)
end
end