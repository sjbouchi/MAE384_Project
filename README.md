clear
clc
h= 1; %the time step is one day
T=100; %number of days
S0= 990; %amount of susceptible people intially
I0= 10; %amount of infected people intially
R0= 0; %amount of recovered people intially
N=S0+I0+R0; %constant, S+R+I
Beta=[.3 1 2];
Gamma=[.1 .1 .2];
%%influenza

[S_influenza,I_influenza,R_influenza]=solve_SIR(Beta(1), Gamma(1), T, h, S0, I0, R0, N);

%%Covid

[S_Covid,I_Covid,R_Covid]=solve_SIR(Beta(2), Gamma(2), T, h, S0, I0, R0, N);

%Measles

[S_Measles,I_Measles,R_Measles]=solve_SIR(Beta(3), Gamma(3), T, h, S0, I0, R0, N);
%generate plots
t = 1:h:T; %time matrix in days

figure(1); %influenza plot
plot(t, S_influenza, 'b',t,I_influenza,'r',t,R_influenza,'g'); %S for influenza
hold on
title('Influenza plot');
xlabel('time in days');
ylabel('population');
legend({'Suseptible population','Infected population','Recovered population'});
hold off

figure(2); %Covid plot
plot(t, S_Covid,'b',t, I_Covid,'r',t,R_Covid,'g'); %S for covid
hold on
title('Covid plot');
xlabel('time in days');
ylabel('population');
legend({'Suseptible population','Infected population','Recovered population'});
hold off

figure(3); %measels plot
plot(t, S_Measles,'b',t,I_Measles,'r',t,R_Measles,'g'); %s for measels
hold on
title('Measels plot');
xlabel('time in days');
ylabel('population');
legend({'Suseptible population','Infected population','Recovered population'});
hold off
The first two plots make sense in the context of the problem as you can see direct correlation between the beta and gamma and the solve functions in the graphs
%) Part 2
clear;
clc;

h = 1; % Time step for fine grid
T = 100; % Number of days
S0 = 990; % Initial susceptible population
I0 = 10; % Initial infected population
R0 = 0; % Initial recovered population
N = S0 + I0 + R0; % Total population
Beta = [.3, 1, 2];
Gamma = [.1, .1, .2];

%% Influenza
[S_h1, I_h1, R_h1] = solve_SIR(Beta(1), Gamma(1), T, h, S0, I0, R0, N);

h_coarse = 2;
[S_h2, I_h2, R_h2] = solve_SIR(Beta(1), Gamma(1), T, h_coarse, S0, I0, R0, N);

% Time vectors
t_h1 = 0:h:T; % h1
t_h2 = 0:h_coarse:T; % h2

% odd days
odd_days = t_h1(mod(t_h1, 2) == 1); 

S_true = zeros(size(odd_days));
I_true = zeros(size(odd_days));
R_true = zeros(size(odd_days));

S_linear_interp = zeros(size(odd_days));
I_linear_interp = zeros(size(odd_days));
R_linear_interp = zeros(size(odd_days));

S_quadratic_interp = zeros(size(odd_days));
I_quadratic_interp = zeros(size(odd_days));
R_quadratic_interp = zeros(size(odd_days));

% Loop over odd days for interpolation
for k = 1:length(odd_days)
    S_true(k) = S_h1(odd_days(k) + 1); 
    I_true(k) = I_h1(odd_days(k) + 1);
    R_true(k) = R_h1(odd_days(k) + 1);
    
    
    idx_prev = max(1, floor(odd_days(k) / 2)); % Index for t1
    idx_curr = min(length(t_h2), idx_prev + 1); % Index for t2
    idx_next = min(length(t_h2), idx_prev + 2); % Index for t3
    
    t1 = t_h2(idx_prev); 
    t2 = t_h2(idx_curr); 
    t3 = t_h2(idx_next); 

   
    if idx_curr < length(t_h2) 
        S_linear_interp(k) = S_h2(idx_curr) + ...
                             (S_h2(idx_next) - S_h2(idx_curr)) * ...
                             (odd_days(k) - t2) / (t3 - t2);

        I_linear_interp(k) = I_h2(idx_curr) + ...
                             (I_h2(idx_next) - I_h2(idx_curr)) * ...
                             (odd_days(k) - t2) / (t3 - t2);

        R_linear_interp(k) = R_h2(idx_curr) + ...
                             (R_h2(idx_next) - R_h2(idx_curr)) * ...
                             (odd_days(k) - t2) / (t3 - t2);
    end

    % Quadratic interpolation 
    if idx_curr < length(t_h2) - 1 
        S_quadratic_interp(k) = S_h2(idx_prev) * ...
                                ((odd_days(k) - t2) * (odd_days(k) - t3)) / ...
                                ((t1 - t2) * (t1 - t3)) + ...
                                S_h2(idx_curr) * ...
                                ((odd_days(k) - t1) * (odd_days(k) - t3)) / ...
                                ((t2 - t1) * (t2 - t3)) + ...
                                S_h2(idx_next) * ...
                                ((odd_days(k) - t1) * (odd_days(k) - t2)) / ...
                                ((t3 - t1) * (t3 - t2));
                            
        I_quadratic_interp(k) = I_h2(idx_prev) * ...
                                ((odd_days(k) - t2) * (odd_days(k) - t3)) / ...
                                ((t1 - t2) * (t1 - t3)) + ...
                                I_h2(idx_curr) * ...
                                ((odd_days(k) - t1) * (odd_days(k) - t3)) / ...
                                ((t2 - t1) * (t2 - t3)) + ...
                                I_h2(idx_next) * ...
                                ((odd_days(k) - t1) * (odd_days(k) - t2)) / ...
                                ((t3 - t1) * (t3 - t2));
                            
        R_quadratic_interp(k) = R_h2(idx_prev) * ...
                                ((odd_days(k) - t2) * (odd_days(k) - t3)) / ...
                                ((t1 - t2) * (t1 - t3)) + ...
                                R_h2(idx_curr) * ...
                                ((odd_days(k) - t1) * (odd_days(k) - t3)) / ...
                                ((t2 - t1) * (t2 - t3)) + ...
                                R_h2(idx_next) * ...
                                ((odd_days(k) - t1) * (odd_days(k) - t2)) / ...
                                ((t3 - t1) * (t3 - t2));
    end
end

%error
E_S_linear = sqrt(mean((S_true - S_linear_interp).^2));
E_I_linear = sqrt(mean((I_true - I_linear_interp).^2));
E_R_linear = sqrt(mean((R_true - R_linear_interp).^2));

E_S_quadratic = sqrt(mean((S_true - S_quadratic_interp).^2));
E_I_quadratic = sqrt(mean((I_true - I_quadratic_interp).^2));
E_R_quadratic = sqrt(mean((R_true - R_quadratic_interp).^2));

%create a table to show the error values
Scol= [E_S_linear, E_S_quadratic]; %columb one
Icol= [E_I_linear, E_I_quadratic]; %columb two
Rcol= [E_R_linear, E_R_quadratic]; %columb three
Table= table(Scol', Icol', Rcol', 'VariableNames',{'Susceptible population','Infected population','Recovered population'}, 'RowNames',{'Linear Interpolation error','Quadratic Interpolation error'});
disp(Table)
Although the errors look large here, it is likely due to the nature of the solution of the differential equations, this can also be due to the large numbers of the functions.
%)Part 3)
clear
clc

% Step one - run non liner SRI model to find "true" data

h = 1; % Days
beta = 0.3;
gamma = 0.1;
T = 30; % Total days
t = 1:T; 
S0=990;
I0=10;
R0=0;
N=S0+I0+R0;

[S, I, R] = solve_SIR(beta, gamma, T, h, S0, I0, R0, N);
% Step two - apply linear least squares model
N = 1000; S0 = 990; Gamma = 0.1; % Parameters given
n = 30; % number of data points


a1 = (n*sum(log(I).*t)-(sum(t)*sum(log(I))))/(n*sum(t.^2)-((sum(t))^2)); % Use linear least squares to calculate a1 = k
a0 = sum(log(I))/n - a1*sum(t)/n; % Use linear least squares to calculate a0 = ln(I(0))

I0 = exp(a0);
b = (a1 + Gamma)*N/S0;

% Step three - repeat step two using only 10 days

Iprime = I(1,1:10); % Create new value of I using only first 10 values
tprime = 1:10;
nprime = 10;

a1prime = (nprime*sum(log(Iprime).*tprime)-(sum(tprime)*sum(log(Iprime))))/(nprime*sum(tprime.^2)-((sum(tprime))^2));
a0prime = sum(log(Iprime))/nprime - a1prime*sum(tprime)/nprime;

I0prime = exp(a0prime);
bprime = (a1prime + Gamma)*N/S0;

%create table including true value and both approximations of I(0) and beta
T = table([10;0.3],[I0;b],[I0prime;bprime],'VariableNames',{'True Value','T = 30','T = 10'},'RowNames',{'I(0)','beta'});

disp(T) % Display table

These estimates make sense within the context of the problem as the initial conditions derived using least squares they are within a reasonable distance t0 I0=10 and beta =.3
% Part 4)
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
