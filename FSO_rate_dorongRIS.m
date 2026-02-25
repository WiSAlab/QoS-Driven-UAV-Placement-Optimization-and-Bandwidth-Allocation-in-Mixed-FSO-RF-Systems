clear all;
%%-----------------General parameter----------------------%%
lamda = 1550*10^(-9); % optical wavelength
L_total = 2500; % tổng khoảng cách
sigma_theta = 3*10^(-5);
sigma_beta  = 2*10^(-4);
sigma_n = 10^-7; % AWGN variance
B_i = 1e9;
a = 10*10^-2;   % user receiver aperture
alpha_coeff = 0.434;
Cn2 = 10^(-14);
divergence_angle = 0.5e-3;
P_t_dBm = 5:5:15;
P_t = 10.^(P_t_dBm/10 -3); % W

L_so = 1500; % RIS vị trí cố định
L_ou = L_total - L_so;

ORIS_sizes = [ 0.6, 0.8, 1.0]; % độ rộng vật lý của RIS (m)
N = 3; % số phần tử RIS cố định
R_all = zeros(length(ORIS_sizes), length(P_t_dBm)); % Kết quả

for j = 1:length(ORIS_sizes)
    ORIS_size = ORIS_sizes(j);
    d_side = ORIS_size / N; % kích thước mỗi phần tử

    %%-----------------Channel model--------------------------%%
    h_l = exp(-alpha_coeff * L_total/1000);
    k_wave = 2*pi/lamda;
    z_L = sqrt(k_wave*a^2/L_total);
    sigma_R = 1.23*Cn2*(2*pi/lamda)^(7/6)*L_total^(11/6);
    sigma_s = sqrt(exp(((0.49*sigma_R^2)/(1+0.18*z_L^2+0.56*sigma_R^(12/5))^(7/6)) + ((0.51*sigma_R^2)/(1+0.9*z_L^2+0.62*sigma_R^(5/6)))^(5/6)) -1);
    alpha = 4.24;
    beta = 2.42;

    w_0 = (2*lamda)/(pi*divergence_angle);
    w_L = w_0*sqrt(1+((lamda*L_so)/(pi*w_0^2))^2);
    v = (a/w_L)*sqrt(pi/2);
    w_Leq = sqrt(w_L^2 * ((sqrt(pi)*erf(v))/(2*v*exp(-v^2))));
    sigma_p = sqrt(sigma_theta^2*L_total^2 + 4*(sigma_beta^2)*L_so^2);
    xi = w_Leq / (2*sigma_p);
    A_0 = (erf(v))^2;

    parfor i = 1:length(P_t_dBm)
        P_over =  3*(1/3) * P_t(i) * exp(-d_side^2/(2*w_L^2));  % Phân tán công suất theo kích thước phần tử
        a_func = 2*P_t(i)^2;
        b_func = 2*P_over^2;
        c_func = sigma_n^2;

        h_channel = @(SNR) sqrt(SNR.*c_func./(a_func - SNR.*b_func));
        h_phay = @(SNR) (c_func*a_func./(2*(a_func - b_func.*SNR).^2)) .* sqrt((a_func - b_func.*SNR)./(c_func.*SNR));
        func_PDF_channel = @(SNR) ((alpha*beta*xi^2)/(A_0*h_l*gamma(alpha)*gamma(beta))) ...
            .* MeijerG({[],[xi^2]},{[xi^2-1, alpha-1, beta-1],[]}, alpha*beta.*h_channel(SNR)./(A_0*h_l));
        fun_SNR = @(SNR) h_phay(SNR).*func_PDF_channel(SNR);

        R_all(j,i) = B_i * integral(@(SNR) log2(1 + SNR) .* fun_SNR(SNR), 0, inf);
    end
end

%%-----------------Plot-----------------------------------%%
figure; hold on;
line_styles = {'-bs','-k^','-md','-ro'};
for j = 1:length(ORIS_sizes)
    semilogy(P_t_dBm, R_all(j,:), line_styles{j}, 'LineWidth', 2, ...
        'DisplayName', sprintf('RIS width = %.1f m', ORIS_sizes(j)));
end
xlabel('Transmit Power (dBm)');
ylabel('FSO Rate R_i (bps)');
% title('Ảnh hưởng độ rộng vật lý RIS đến tốc độ dữ liệu');
legend('Location', 'best');
grid on;
