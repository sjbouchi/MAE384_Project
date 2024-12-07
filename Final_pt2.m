%% Part 2: Linear and Quadratic Interpolation on Odd Days
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