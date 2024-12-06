clear
clc


T = 30;
h = 0.1;
N = 1000;
S0 = 990;
I0 = 10;
R0 = 0;

Beta0 = 0.3;
Gamma = 0.1;
A = 5;
omega1 = 2 * pi / 1; 
omega2 = 2 * pi / (365 / 100);

t = 0:h:T;
S = zeros(1, length(t));
I = zeros(1, length(t));
R = zeros(1, length(t));

S(1) = S0;
I(1) = I0;
R(1) = R0;

for i = 2:length(t)
    Beta_t = Beta0 * (1 + A * sin(omega1 * t(i)));
    dS = -(Beta_t / N) * S(i-1) * I(i-1);
    dI = (Beta_t / N) * S(i-1) * I(i-1) - Gamma * I(i-1);
    dR = Gamma * I(i-1);

    S(i) = S(i-1) + dS * h;
    I(i) = I(i-1) + dI * h;
    R(i) = R(i-1) + dR * h;
end

figure;
plot(t, S, 'b', t, I, 'r', t, R, 'g');
title('SIR Model with Daily Periodic Transmission Rate');
xlabel('Time (days)');
ylabel('Population');
legend('S(t)', 'I(t)', 'R(t)');

N_samples = length(I);
f = (0:N_samples/2-1) / (h * N_samples);
I_fft = abs(fft(I));
I_fft_half = I_fft(1:N_samples/2);

figure;
plot(f, I_fft_half);
title('Frequency Spectrum of I(t) with Daily Periodicity');
xlabel('Frequency (1/day)');
ylabel('Magnitude');

S = zeros(1, length(t));
I = zeros(1, length(t));
R = zeros(1, length(t));

S(1) = S0;
I(1) = I0;
R(1) = R0;

for i = 2:length(t)
    Beta_t = Beta0 * (1 + A * sin(omega2 * t(i)));
    dS = -(Beta_t / N) * S(i-1) * I(i-1);
    dI = (Beta_t / N) * S(i-1) * I(i-1) - Gamma * I(i-1);
    dR = Gamma * I(i-1);

    S(i) = S(i-1) + dS * h;
    I(i) = I(i-1) + dI * h;
    R(i) = R(i-1) + dR * h;
end

figure;
plot(t, S, 'b', t, I, 'r', t, R, 'g');
title('SIR Model with ~3-Day Periodic Transmission Rate');
xlabel('Time (days)');
ylabel('Population');
legend('S(t)', 'I(t)', 'R(t)');

% Fourier Transform for ~3-day periodicity
I_fft = abs(fft(I));
I_fft_half = I_fft(1:N_samples/2);

figure;
plot(f, I_fft_half);
title('Frequency Spectrum of I(t) with ~3-Day Periodicity');
xlabel('Frequency (1/day)');
ylabel('Magnitude');
