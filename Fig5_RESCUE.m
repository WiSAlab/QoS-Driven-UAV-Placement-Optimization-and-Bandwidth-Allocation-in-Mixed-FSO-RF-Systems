clear;

% General Parameters
% FSO transmission power (~ W)
P_t     = 0.1; 

% Divergence angle (~ rad)
theta_g = 1e-3;

% Receiver radius (~ m)
r       = 0.05;

% Receiver sensitivity (~ photons/bits)
 N_b     = 67885;
%  N_b = 100000;

% FSO wavelength (~m)
lambda  = 1550e-9;

% Visible distance (~ m)
v       = 5e3;

% Disaster area radius (~ m)
area_radius = 2000;

% Min/Max altitude of UAV (~ m)
h_max   = 1000;
h_min   = 100;

% Carrier frequency (~ Hz)
f_c     = 2e9;

% Environment index 
b_idx    = 9.61;
beta_idx = 0.16;

% Average excessive pathloss in LoS and NLoS (~ dB)
xi_LoS  = 1;
xi_NLoS = 20;
% xi_LoS  = 10^(1/10);
% xi_NLoS = 10^(20/10);

% Noise power spectral density (~ dBm/Hz)
N_0     = -104;

% UAV downlink transmission power (~ dBm)
P_t_down    = 20;

% Available bandwidth for each UAV (~ Hz)
B_i     = 1e9;

% Pathloss requirement of MUs (~ dB)
eta_th  = 96;

% UAV & BS - Initialization
num_UAVs   = 4;    % number of UAVs
UAV_pos = [[500, 500]; [1500, 500]; [1500, 1500]; [500, 1500]]; 
% BS_pos = [[-500, 1000]; [1000, -500]; [2500, 1000]; [1000, 2500]];
BS_pos =[2500,1000]; %OIRS position

UAV = struct;
for i = 1:num_UAVs
    UAV(i).id = i;
    UAV(i).pos = UAV_pos(i, :); % position of the UAV
%     UAV(i).BS = BS_pos(i, :); % position of the BS associated with the UAV
    UAV(i).BS = BS_pos;
    UAV(i).h = 0;   % altitude of the UAV
    UAV(i).d = 0;   % euclid distance from the BS to the UAV
    UAV(i).B = B_i; % total bandwidth allocated to the UAV
    UAV(i).B_access = 0; % total bandwidth offered to MUs by the UAV
    UAV(i).rate_fso = 0; % backhaul rate of the UAV
    UAV(i).rate_access = 0; % total rate offered to MUs by the UAV
    UAV(i).MU_id = [];  % ID of the MUs supported by the UAV
end


% Mobile Users - Initialization
MU_density_list = 20:30; % Các mật độ cần mô phỏng
results = zeros(length(MU_density_list), 1); % Lưu tỷ lệ MU được phục vụ
for density_idx = 1:length(MU_density_list)
    MU_density = MU_density_list(density_idx);
    
% % MU_density = 0.5; % average MUs per location
% MU_pos = zeros(10*10*1, 2);
% num_MUs = 0; 
% 
% MU_center = zeros(10*10, 2);
% x_axis = -10;
% for i = 1:10
%     x_axis = x_axis + 200;
%     y_axis = -10;
%     for j = 1:10
%         y_axis = y_axis + 200;
%         MU_center(10*(i-1) + j, :) = [x_axis, y_axis];  
%         [xy, numPoints] = MU_position_gen(10, MU_density, [x_axis, y_axis]);
%         MU_pos(num_MUs+1 : num_MUs + numPoints, :) = xy;
%         num_MUs = num_MUs + numPoints;
%     end
% end

% MU_pos(num_MUs+1 : end, :) = [];
% Khởi tạo MU 
    MU_center = zeros(10*10*1, 2);
    num_MUs = 0;
    x_axis = -10;
for i = 1:10
    x_axis = x_axis + 200;
    y_axis = -10;
    for j = 1:10
        y_axis = y_axis + 200;
        MU_center(10*(i-1) + j, :) = [x_axis, y_axis];  
        [xy, numPoints] = MU_position_gen(10, MU_density, [x_axis, y_axis]);
        MU_pos(num_MUs+1 : num_MUs + numPoints, :) = xy;
        num_MUs = num_MUs + numPoints;
    end
end
% MU_rate = normrnd(3e6, 1e6, [1, num_MUs]);
MU_rate = rand(1, num_MUs) .* 20e6 + 1e6; 

MU = struct;
for i = 1:num_MUs
    MU(i).id = i;
    MU(i).pos = MU_pos(i,:);
    MU(i).rate = MU_rate(i);
    MU(i).pathloss = 0;
    MU(i).bandwidth = 0;
    MU(i).UAV = 0;
end

% Main function
r_opt = zeros(1, num_UAVs);

a     = zeros(num_UAVs, num_MUs);

for UAV_idx = 1:num_UAVs
    
%     Initial optimal altitude of UAV
    func1 = @(theta) pi/(9*log(10)) .* tan(theta) + b_idx*beta_idx*(xi_LoS - xi_NLoS)...
        .* exp(-beta_idx.*(180/pi*theta - b_idx)) ./ (b_idx.*exp(-beta_idx.*(180/pi*theta - b_idx)) +1).^2;
    
    theta_opt = bisection(func1, 0, pi/2, 0); % in rad
    
    func2 = @(r) (xi_LoS - xi_NLoS) ./ (1 + b_idx.*exp(-beta_idx.*(180/pi*theta_opt - b_idx))) ...
        + 20*log10(r/cos(theta_opt) ) + 20*log10(4*pi*f_c/3e8) + xi_NLoS;
    
    r_opt(UAV_idx) = bisection(func2, 0, 2000, eta_th);
    
    UAV(UAV_idx).h = r_opt(UAV_idx) * tan(theta_opt);    

%     UAV: Altitude adjustment
    delta_min = 10;
    delta_max = 100;
    delta = delta_min;

    n_MU_supported = zeros(1,3);    % [h_cur, h_plus, h_minus]

    while delta <= delta_max
        h_minus = max(UAV(UAV_idx).h - delta, h_min);
        h_plus = min(UAV(UAV_idx).h + delta, h_max);

        h_arr = [UAV(UAV_idx).h, h_plus, h_minus];

        for h_idx = 1:length(h_arr)
            UAV(UAV_idx).h = h_arr(h_idx);
%             Backhaul rate
        UAV(UAV_idx).d = euclid_dist(UAV(UAV_idx).pos, UAV(UAV_idx).BS, 0, UAV(UAV_idx).h);
        sigma = atmos_attenuation(lambda, v);
        UAV(UAV_idx).rate_fso = rate_fso(P_t, lambda, N_b, r, theta_g, UAV(UAV_idx).d, sigma);
    
%             Access Link: Pathloss
        for MU_idx = 1:num_MUs
            MU(MU_idx).pathloss = pathloss_user(MU(MU_idx).pos, UAV(UAV_idx).pos, UAV(UAV_idx).h, ...
                    f_c, xi_LoS, xi_NLoS, b_idx, beta_idx);
        end
    
%             Access Link: required bandwidth
        MU_covered = MU([MU.pathloss] <= eta_th);

        if isempty(MU_covered)
            n_MU_supported(h_idx) = 0;
            continue;
        end
        
        % MU_bandwidth  = MU_covered_rate ./ log2(1 + 10^(P_t_down/10-3)/10^(N_0/10-3) .* 10.^(-MU_covered_pathloss ./ 10));
        MU_bandwidth = [MU_covered.rate] ./ log2(1 + 10^(P_t_down/10-3)/10^(-20.4) .* 10.^(-[MU_covered.pathloss] ./ 10));
       
        for i = 1:length(MU_covered)
            MU_covered(i).bandwidth = MU_bandwidth(i);
        end
%             Bandwidth allocation and MU association
        [~, Index] = sort([MU_covered.rate]);
        MU_sorted = MU_covered(Index);
        
        i = 0;
        MU_idx_arr = zeros(1, length(MU_covered));
        UAV(UAV_idx).B_access = 0;
        UAV(UAV_idx).rate_access = 0;

        while 1
            i = i+1;
            
            % if the current MU is not associated with any UAVs
%             if MU(MU_sorted(i).id).UAV == 0
            if i > length(MU_sorted)
                 break; % Ngăn lỗi truy cập ngoài phạm vi
            end

            if MU(MU_sorted(i).id).UAV == 0
                new_B_access = UAV(UAV_idx).B_access + MU_sorted(i).bandwidth;
                new_rate_access = UAV(UAV_idx).rate_access + MU_sorted(i).rate;
            end
    
            if (i > length(MU_sorted)) || (new_B_access > UAV(UAV_idx).B) ...
                || (new_rate_access > UAV(UAV_idx).rate_fso) 
                break;
            end
            
            MU_idx_arr(i) = MU_sorted(i).id;
            UAV(UAV_idx).B_access = new_B_access;
            UAV(UAV_idx).rate_access = new_rate_access;
            
        end
%         MU_idx_arr(i:end) = []; 
        if ~isempty(MU_idx_arr)
            MU_idx_arr(i:end) = [];   
        else
            MU_idx_arr = []; % Đảm bảo không bị lỗi khi rỗng
        end
        
        % Assign the IDs of the highest number of MUs to the UAV
        if length(MU_idx_arr) >= length(UAV(UAV_idx).MU_id)
            UAV(UAV_idx).MU_id = MU_idx_arr;
        end
        
%         Altitude Adjustment
        n_MU_supported(h_idx) = i-1;
        end 
        
        % Update the alitude of the UAV
        if max(n_MU_supported) == n_MU_supported(1)
            UAV(UAV_idx).h = h_arr(1);
            delta = delta + delta_min;
        elseif max(n_MU_supported) == n_MU_supported(2)
            UAV(UAV_idx).h = h_arr(2);
            delta = delta_min;
        elseif max(n_MU_supported) == n_MU_supported(3)
            UAV(UAV_idx).h = h_arr(3);
            delta = delta_min; 
        end

    end % end updating the altitude loop
    
    % Assign UAV_id that associate to the MU
    for i = 1:length(UAV(UAV_idx).MU_id)
        id = UAV(UAV_idx).MU_id(i);
        MU(id).UAV = UAV_idx;
    end

end


num_MU_served = sum([MU.UAV] > 0);
results(density_idx) = num_MU_served / num_MUs * 100;
end

figure;
plot(MU_density_list, results, 'b-', 'LineWidth', 2);
xlabel('Average number of MUs per 20m x 20m location', 'FontSize', 12);
ylabel('Fraction of served MU(%)', 'FontSize', 12);
title('Fraction of served MUs data rate requirements', 'FontSize', 14);
grid on;
xticks(MU_density_list);
xlim([min(MU_density_list), max(MU_density_list)]);
ylim([0 100]);
