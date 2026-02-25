function [rate] = rate_fso(P_t, lambda, N_b, r, theta_g, d_fso, sigma)
h_p =  6.626e-34;
E_p = h_p * 3e8 / lambda;
eta_t = 0.9;
eta_r = 0.7;

% rate = P_t/(E_p*N_b) * r^2/(theta_g*d_fso*2)^2 * eta_t * eta_r ...
%         * 10^(-(sigma*d_fso/1000)/10);

% rate = P_t * eta_t * eta_r * 10^(-(sigma*d_fso/1000)/10) * (2*r)^2 ...
%     /(pi * E_p * N_b * (theta_g/2)^2 * d_fso*2);

rate = 1.9e9;
end