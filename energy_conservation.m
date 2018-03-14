% finding the maximum amplitude of wave two using
% energy conservation, using wave one

global wave_number; % calling wave_number from dispersion script
global sigma_sq; % calling sigma_sq from dispersion script
global h;
global g;

sigma = sqrt(sigma_sq);
a_1 = 0.5; % amplitude of deep water wave [m]

C_p1 = sigma/wave_number; % group velocity
n = 0.5; % because in deep water; in shallow n = 1
C_g1 = C_p1*n;
C_g2 = sqrt(9.81*h);

amp_2_var = a_1*(sqrt(C_g1/C_g2))
amp_2 = vpa(amp_2_var)

u_velocity = amp_2*sigma/(wave_number*h)


