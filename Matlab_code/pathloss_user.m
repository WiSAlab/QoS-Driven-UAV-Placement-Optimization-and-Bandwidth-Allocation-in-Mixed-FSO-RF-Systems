function [eta_avr] = pathloss_user(MU_pos, UAV_pos, h_UAV, f_c, xi_LoS, xi_NLoS, ...
                    b_idx, beta_idx)

% d = euclid_dist(MU_pos(1,:), UAV_pos(UAV_idx, :), 0, h_UAV(UAV_idx));
% eta_LoS = 20*log10(4*pi*f_c*d/3e8) + xi_LoS; 
% eta_NLoS = 20*log10(4*pi*f_c*d/3e8) + xi_NLoS;
% p_LoS = 1 / (1 + b_idx*exp(-beta_idx*(180/pi*asin(h_UAV(UAV_idx) /d ) - b_idx)) );
% eta_avr = p_LoS * eta_LoS + (1 - p_LoS)*eta_NLoS;

d = euclid_dist(MU_pos, UAV_pos, 0, h_UAV);
eta_LoS = 20*log10(4*pi*f_c*d/3e8) + xi_LoS; 
eta_NLoS = 20*log10(4*pi*f_c*d/3e8) + xi_NLoS;
p_LoS = 1 / (1 + b_idx*exp(-beta_idx*(180/pi*asin(h_UAV /d ) - b_idx)) );
eta_avr = p_LoS * eta_LoS + (1 - p_LoS)*eta_NLoS;

end

