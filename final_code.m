clear
clc
% part one: ode solver
h= 1; %the time step is one day
T=100; %number of days
%influenza
S= zeros(1,T); %Intialize susceptible individuals matrix
I= zeros(1,T); %Intialize infected individuals matrix
R= zeros(1,T); %intialize recovered individuals matrix
S(1)= 990; %amount of susceptible people intially
I(1)= 10; %amount of infected people intially
R(1)= 0; %amount of recovered people intially
N=1000; %constant, S+R+I

Beta= .3; %transmission rate
Gamma= .1; %recovery rate
dS= @(I,S) -(Beta/N)*S*I; %ode for susceptible
dI= @(S,I) (Beta/N)*S*I-(Gamma*I); %ode for infected rate
dR= @(I) Gamma*I; %ode for recovery rate
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
%realocate values to different matrix for influenza
S_influenza=S;
I_influenza=I;
R_influenza=R;

%Covid
S= zeros(1,T); %Intialize susceptible individuals matrix
I= zeros(1,T); %Intialize infected individuals matrix
R= zeros(1,T); %intialize recovered individuals matrix
S(1)= 990; %amount of susceptible people intially
I(1)= 10; %amount of infected people intially
R(1)= 0; %amount of recovered people intially
N=1000; %constant, S+R+I

Beta= 1; %transmission rate
Gamma= .1; %recovery rate
dS= @(I,S) -(Beta/N)*S*I; %ode for susceptible
dI= @(S,I) (Beta/N)*S*I-(Gamma*I); %ode for infected rate
dR= @(I) Gamma*I; %ode for recovery rate
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
%realocate values to different matrix for influenza
S_Covid=S;
I_Covid=I;
R_Covid=R;

%Measles
S= zeros(1,T); %Intialize susceptible individuals matrix
I= zeros(1,T); %Intialize infected individuals matrix
R= zeros(1,T); %intialize recovered individuals matrix
S(1)= 990; %amount of susceptible people intially
I(1)= 10; %amount of infected people intially
R(1)= 0; %amount of recovered people intially
N=1000; %constant, S+R+I

Beta= 2; %transmission rate
Gamma= .2; %recovery rate
dS= @(I,S) -(Beta/N)*S*I; %ode for susceptible
dI= @(S,I) (Beta/N)*S*I-(Gamma*I); %ode for infected rate
dR= @(I) Gamma*I; %ode for recovery rate
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
%realocate values to different matrix for influenza
S_Measels=S;
I_Measels=I;
R_Measels=R;

%generate plots
t= (1:h:T); %time matrix in days
figure(1) %influenza plot
plot(t, S_influenza, 'b') %S for influenza
hold on
plot(t,I_influenza,'r') %I for influenza
plot(t,R_influenza,'g')%R for influenza
title('Influenza plot')
xlabel('time in days')
ylabel('population')
legend({'Suseptible population','Infected population','recovered population'})
hold off

figure(2) %Covid plot
plot(t, S_Covid,'b') %S for covid
hold on
plot(t, I_Covid,'r')%I for covid
plot(t,R_Covid,'g')%R for covid
title('Covid plot')
xlabel('time in days')
ylabel('population')
legend({'Suseptible population','Infected population','recovered population'})
hold off

figure(3) %measels plot
plot(t, S_Measels,'b') %s for measels
hold on
plot(t,I_Measels,'r')%I for meassels
plot(t,R_Measels,'g') %R for meassels
title('Measels plot')
xlabel('time in days')
ylabel('population')
legend({'Suseptible population','Infected population','recovered population'})
hold off