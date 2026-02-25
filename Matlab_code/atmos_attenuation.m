function sigma_dB = atmos_attenuation(lambda, v)

V = v/1000;
if V > 50
    q = 1.6;
elseif V > 6
    q = 1.3;
elseif V > 1
    q = 0.16*V + 0.34;
elseif V > 0.5
    q = V - 0.5;
else 
    q = 0;
end

sigma_dB = (3.912/V)*(lambda*10^9/550)^(-q);  %attenuation coefficent (dB/km) 
% sigma    = sigma_dB/(10^4*log10(exp(1)));       %Unit: m^-1
% h_l       = exp(-sigma*H_cl*(1/cos(xi)));        %Could attenuation Calculation 
end
