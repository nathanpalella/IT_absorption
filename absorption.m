%% Tube geometry and preliminary variable configuration
x1 = 1;     % in, distance between reference mic and beginning of specimen
s = 4;      % in, distance between mics 
d = 2;      % in, diameter of tube

[mic_a,Fs] = audioread('mic_a.wav');  % in the standard position
mic_b = audioread('mic_b.wav');
L = length(mic_a); 

FFT_a = fft(mic_a);
FFT_b = fft(mic_b);
S11 = real(FFT_a).*FFT_a;
S12 = real(FFT_b).*FFT_a;

H12 = S12./S11;

d_meter = .0254*d;    % conversion to IS units
s_meter = .0254*s;
x1_meter = .0254*x1;

T = 313;    % k, ambient air temperature
T_0 = 293;  % K, standard air pressure
p = 97;     % kPa, ambient air pressure in kilopascals
p_0 = 101.325;   % kPa, standard air pressure
rho_0 = 1.186;   % kg/m^3, standard air density
rho = rho_0 * ( (p*T_0) / (p_0*T)); % kg/m^3, ambient air density

c0 = 343.2 * sqrt(T/T_0); % m/s, calculated speed of sound
% Maximum and minimum frequencies
fu = (0.58*c0)/d_meter;    % Hz, maximum frequency
fl = 100;  % Hz, chosen based on recording equipment

% Reminder
% k0 = 2*pi*f_a/c0; % wave numbers
% y0 = c0./f_a;      % m, wavelengths

%% Calibration using a predetermined calibration sample

H12_1 = FFT_A ./ FFT_B;

mic_a_swap = audioread('mic_a_swap.wav');  % in the switched position
mic_b_swap = audioread('mic_a_swap.wav');
FFT_A_swap = fft(mic_a_swap); FFT_B_swap = fft(mic_b_swap);

H12_2 = (FFT_A_swap ./ FFT_B_swap);

Hc = (H12_1 ./ H12_2) .^ 0.5;     % correction factor with specimen in situ

H12_c = H12./Hc;  % Correct the original mic placement test

HI = zeros(L/2,1);  % preallocate memory for efficiency
HR = zeros(L/2,1);
k0 = zeros(L/2,1);
for x = 1:1:L/2
    k0(x,1) = ((2*pi*Fs)/c0);
end
for y = 1:1:L/2
    HI(y,1) = exp(-i*k0(y,1)*s_meter);
    HR(y,1) = exp(i*k0(y,1)*s_meter);
end

r = ( (H12_c - HI) ./ (HR - H12_c) ) .* exp(2.*i.*k0.*x1_meter);
a = 1 - abs(r.^2);
a = abs(a);

f = 1:Fs/L:Fs/2;
%plot results 
freq = zeros(L/2,1); 
for z=1:1:L/2
    freq(z,1) = z*Fs-1; 
end

plot(freq,a)