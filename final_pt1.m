clear
clc
% part one: ode solver
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
t = 0:h:T; %time matrix in days
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
