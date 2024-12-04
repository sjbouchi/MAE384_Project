clear
clc

% Step one - run non liner SRI model to find "true" data

h = 1; % Days
beta = 0.3;
gamma = 0.1;
T = 30; % Total days
t = 1:T; 

S = zeros(1,T); % Preallocate Vector for S, I, R
I = zeros(1,T);
R = zeros(1,T);
S(1) = 990; % initial susceptible, infected, and recovered population 
I(1) = 10;
R(1) = 1;

N = S(1) + I(1) + R(1);

dS= @(I,S) -(beta/N)*S*I; %ode for susceptible
dI= @(S,I) (beta/N)*S*I-(gamma*I); %ode for infected rate
dR= @(I) gamma*I; %ode for recovery rate
k=2; %intialize a counter

while k<= T
    %start with k values for S
    k1_S= dS(I(k-1),S(k-1)); 
    k2_S= dS(I(k-1)+.5*h,S(k-1)+.5*k1_S*h);
    k3_S= dS(I(k-1)+.5*h,S(k-1)+.5*k2_S*h);
    k4_S= dS(I(k-1)+h,S(k-1)+k3_S*h);
    S(k)= S(k-1) + (1/6)*(k1_S+ 2*k2_S + 2*k3_S + k4_S)*h; %values for S
    %Now do k values for I
    k1_I= dI(S(k-1),I(k-1)); 
    k2_I= dI(S(k-1)+.5*h,I(k-1)+.5*k1_I*h);
    k3_I= dI(S(k-1)+.5*h,I(k-1)+.5*k2_I*h);
    k4_I= dI(S(k-1)+h,I(k-1)+k3_I*h);
    I(k)= I(k-1) + (1/6)*(k1_I+ 2*k2_I + 2*k3_I + k4_I)*h; %values for I
    %now do k values for R
    k1_R= dR(I(k-1));
    k2_R= dR(I(k-1)+.5*h);
    k3_R= dR(I(k-1)+.5*h);
    k4_R= dR(I(k-1)+h);
    R(k)= R(k-1) + (1/6)*(k1_R + 2*k2_R + 2*k3_R + k4_R); %values for R
    k=k+1; %update counter
end

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
