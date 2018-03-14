T = 12; %s
global sigma_sq;
sigma_sq = ((2*pi)/(T))^2;

global g;
g = 9.81; %m/s^2
global h;
h = 2; %m

% solving dispersion equation for k
% and returning a positive value
% sigma^2 = g*k*tanh(k*h)

syms k
eqn = g*k*tanh(h*k) == sigma_sq;
global wave_number;
wave_number = solve(eqn, k);

if wave_number < 0
    wave_number = wave_number * (-1)
else
    wave_number = wave_number
end

% calculating wavelength
lam = (2*pi)/(wave_number);
lambda = vpa(lam)

% determining deep or shallow water wave
cond = h/lambda
if cond < 0.05
    disp('This is a shallow water wave.')
elseif cond > 0.5
    disp('This is a deep water wave.')
else 
    disp('This is a transitional wave.')
end