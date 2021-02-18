%% Create the supply of Cell N+1
% this is done by create a complementary of the density that tells us when
% there should be a drop in the  Supply

%% Interpolate the density to find a qualitative description
% interpolation with piece-wise linear function
% divide the data in macro sections of length interval
interval = 25;%200 % the bigger the less accurate the interpolation
extremum_points = 1:interval:length(density);
% add the last index of the vector density if it not already in
if extremum_points(end) ~= length(density)
    extremum_points = [extremum_points length(density)];
end
avg_density = zeros(1,length(extremum_points)+1);
index_inter_middle  = zeros(1,length(extremum_points)+1); % index in the vector of density
index_inter_middle(1) = 1; index_inter_middle(end) = extremum_points(end);
for j = 1 : length(extremum_points)-1
    % find the middle index of the inteval "j"
    index_inter_middle(j+1) = extremum_points(j)+round((extremum_points(j+1)-extremum_points(j))/2);
    avg_density(j+1) = mean(density(extremum_points(j):extremum_points(j+1)));
end
% starting and end values doubled 
avg_density(1) = avg_density(2); avg_density(end) = avg_density(end-1); 
%Interpolate the values between two points in "avg_density" in order to
%make it for every time stamps
pwl_density = zeros(size(density)); %(piece-wise-linear  density)
for j = 1 : length(index_inter_middle)-1
    index_pre_interp = index_inter_middle(j:j+1);
    index_post_interp = index_inter_middle(j):index_inter_middle(j+1)-1;     
    pwl_density(index_post_interp) = interp1(index_pre_interp,avg_density(j:j+1),index_post_interp);
end
pwl_density(end) = avg_density(end);

%% Create the supply funciton
% flip the data 
tmp1 = -pwl_density;
% desired amplitude in percentage wrt to the maximum cell capacity
supply_amplitude = 0.9;%0.60;
% compute the scaling factor
tmp1_scaling = supply_amplitude*CTM_param.q_max(n)/(max(tmp1)-min(tmp1));
% correct the amplitude
tmp2 = tmp1_scaling*tmp1;
% shift up to have the maximum around the value of q_max^{N}
supply_N_plus = tmp2+max(tmp2)+1.1*CTM_param.q_max(n);

supply_N_plus(5000:end) = 0.82*CTM_param.q_max(n)/(max(tmp1)-min(tmp1))*tmp1(5000:end)+max(tmp2)+1.1*CTM_param.q_max(n);




