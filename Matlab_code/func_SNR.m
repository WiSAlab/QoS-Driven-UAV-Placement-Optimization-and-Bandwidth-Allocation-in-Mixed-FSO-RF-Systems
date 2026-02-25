% This program calculates the Outage probabbility of OIRS-assisted
% multiuser FSO communications
clear all;
%%-----------------General parameter----------------------%%
lamda = 1550*10^(-9); % the optical wavelength
L_total = 2500; % The total distance
% L_so = 0.5*L_total; % the distance from source to OIRS
% L_ou = 0.5*L_total; % the distance from OIRS to destination
L_so = 1500;
L_ou = 1000;
sigma_theta = 3*10^(-5);%3*10^(-5);
sigma_beta  = 2*10^(-4);
sigma_n = 10^-7; % the AWGN variance
ORIS_size = 1; % the OIRS size
N = 3; % the number of adjacent element
d_side = (1/2)*ORIS_size ; % the side length of one OIRS element
% r_p = 10*10^-2;
a = 10*10^-2;   % the user radius
B_i = 1e9;

%%-----------------Channel model--------------------------%%

% Propagation loss
alpha_coeff = 0.434;
 h_l = exp(-alpha_coeff * L_total/1000);
%  h_l = 1;

% Atmospheric turbulence
Cn2 = 10^(-14);
sigma_R = 1.23*Cn2*(2*pi/lamda)^(7/6)*L_total^(11/6); %1.23*(2*pi/lambda)^(7/6)*Cn2*Lfso^(11/6);
k_wave = 2*pi/lamda;
z_L = sqrt(k_wave*a^2/L_total);
sigma_s = sqrt(exp(((0.49*sigma_R^2)/(1+0.18*z_L^2+0.56*sigma_R^(12/5))^(7/6)) + ((0.51*sigma_R^2)/(1+0.9*z_L^2+0.62*sigma_R^(5/6)))^(5/6)) -1);
alpha = 4.24;%(exp((0.49.*sigma_R.^2)./((1+1.11.*sigma_R.^(12/5)).^(7/6)))-1).^(-1);
beta = 2.42;%(exp((0.51.*sigma_R.^2)./((1+0.69.*sigma_R.^(12/5)).^(5/6)))-1).^(-1);

% Pointing error

divergence_angle = 0.5*10^-3; % the divergence angle
% divergence_angle = source_divergence_angle*(1 + L_so/L_ou);
w_0 = (2*lamda)/(pi*divergence_angle); % the beam-waist at the source
w_L = w_0*sqrt(1+((lamda*L_so)/(pi*w_0^2))^2); % the beam-waist at the distance of L
v = (a/w_L)*sqrt(pi/2);
w_Leq = sqrt(w_L^2 .* ((sqrt(pi)*erf(v))/(2*v*exp(-v^2))));
sigma_p = sqrt(sigma_theta^2*L_total^2 + 4*(sigma_beta^2)*L_so^2);
xi = w_Leq / (2*sigma_p);
A_0 = (erf(v))^2; % the fraction of the collected power

%%-----------------Outage probability---------------------%%

P_t_dBm = 5:5:25;    % the transmit power in dBm
P_t = 10.^(P_t_dBm/10 -3); % Converting from dBm to watt
% N_interfer
R_i = zeros(1, length(P_t_dBm));
parfor i=1:length(P_t_dBm)
    P_over =  3*(1/3) * P_t(i) * exp(-d_side^2/(2*w_L^2)); % Interference optical power
    % P_over_dBm(i) = 10*log10(P_over) + 30;
% Using Log-normal
    % P1 = ((xi^2*P_t(i)^2*sigma_n^2)/(2*(A_0*h_l)^(xi^2)));
    % P2 = @(SNR) ((SNR.*sigma_n^2).^((xi^2)/2 -1).*exp(2*sigma_s^2*xi^2*(1+xi^2)))./((2*P_t(i)^2 - 2*SNR*P_over^2).^((xi^2)/2+1)); %
    % P3 = @(SNR) erfc(real((log((1/(A_0*h_l)).*sqrt((SNR.*sigma_n^2)./(2*P_t(i)^2 - 2.*SNR.*P_over^2))) + 2*sigma_s^2*xi^2*(1+xi^2))/sqrt(8*sigma_s^2)));
    % outage_result(i) = integral(@(SNR) P1.*P2(SNR).*P3(SNR), 0, 5);
% Using Gamma-Gamma: MeijerG({[a_1,...a_n],[a_n+1,...a_p]},{[b_1,...b_m],[b_m+1,...b_q]},z)
    a_func = 2*P_t(i)^2;
    b_func = 2*P_over^2;
    c_func = sigma_n^2;
    % h = @(SNR) sqrt((SNR.*c_func^2)./(a_func + SNR.*b_func));
    % P1 = @(SNR) (c_func*a_func)./((2*(a_func-b_func.*SNR).^2).*sqrt((c_func.*SNR)./(a_func-b_func.*SNR)));
    % SNR_ave = ((P_t(i))^2)/( (P_over)^2 + sigma_n^2);
    % func_SNR = @(SNR) ((xi^2)./(2.*SNR*gamma(alpha)*gamma(beta)))...
    %     .*MeijerG({[],[xi^2+1]},{[xi^2, alpha, beta],[]}, (alpha*beta.*xi^2.*sqrt(SNR))/(sqrt(SNR_ave)*(xi^2+1)));
    % outage_result(i) = integral(func_SNR, 0, 5);
    h_channel = @(SNR) sqrt(SNR.*c_func./(a_func-SNR.*b_func));
    h_phay =@(SNR) (c_func*a_func./(2*(a_func-b_func.*SNR).^2)).*sqrt((a_func - b_func.*SNR)/(c_func.*SNR));
    func_PDF_channel = @(SNR) ((alpha*beta*xi^2)/(A_0*h_l*gamma(alpha)*gamma(beta)))...
        .*MeijerG({[],[xi^2]},{[xi^2-1, alpha-1, beta-1],[]}, alpha*beta.*h_channel(SNR)./(A_0*h_l));
    fun_SNR = @(SNR) h_phay(SNR).*func_PDF_channel(SNR);
%     outage_result(i) = integral(fun_SNR, 0, 10^(0/10));
%       integrand = @(SNR) log2(1 + SNR) .* fun_SNR(SNR);
    R_i(i) = B_i * integral(@(SNR) log2(1 + SNR) .* fun_SNR(SNR), 0, inf);

end

%%-----------------Plot-----------------------------------%%
figure;
semilogy(P_t_dBm, R_i,'-ro');%,P_t, P_t,'r-o');
xlabel('Transmit power (dBm)');
ylabel('Rate FSO');
title('Ảnh hưởng của kênh truyền FSO tới rate FSO R_i');