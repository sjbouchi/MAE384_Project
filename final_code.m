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

%%Part 2 interpolation
%run part one but with larger step size
h= 2; %the time step is one day
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
k=4; %intialize a counter
S(2)=S_influenza(2);
I(2)= I_influenza(2);
R(2)= R_influenza(2); %Intialize 
while k<= T
    %start with k values for S
    k1_S= dS(I(k-2),S(k-2)); 
    k2_S= dS(I(k-2)+.5*h,S(k-2)+.5*k1_S*h);
    k3_S= dS(I(k-2)+.5*h,S(k-2)+.5*k2_S*h);
    k4_S= dS(I(k-2)+h,S(k-2)+k3_S*h);
    S(k)= S(k-2) + (1/6)*(k1_S+ 2*k2_S + 2*k3_S + k4_S)*h; %values for S
    %Now do k values for I
    k1_I= dI(S(k-2),I(k-2)); 
    k2_I= dI(S(k-2)+.5*h,I(k-2)+.5*k1_I*h);
    k3_I= dI(S(k-2)+.5*h,I(k-2)+.5*k2_I*h);
    k4_I= dI(S(k-2)+h,I(k-2)+k3_I*h);
    I(k)= I(k-2) + (1/6)*(k1_I+ 2*k2_I + 2*k3_I + k4_I)*h; %values for I
    %now do k values for R
    k1_R= dR(I(k-2));
    k2_R= dR(I(k-2)+.5*h);
    k3_R= dR(I(k-2)+.5*h);
    k4_R= dR(I(k-2)+h);
    R(k)= R(k-2) + (1/6)*(k1_R + 2*k2_R + 2*k3_R + k4_R); %values for R
    k=k+2; %update counter
end
%interpolation
Si=S; %I am going to use the Matrix Si for linear interpolation
Ii=I;
Ri=R;
c=3; %counter for interpolation
while c<= 99
Si(c)=(((c-(c+1))/((c-1)-(c+1)))*Si(c-1)) + (((c-(c-1))/((c+1)-(c-1)))*Si(c+1)); %linear interpolation lagrange form
Ii(c)=(((c-(c+1))/((c-1)-(c+1)))*Ii(c-1)) + (((c-(c-1))/((c+1)-(c-1)))*Ii(c+1));
Ri(c)=(((c-(c+1))/((c-1)-(c+1)))*Ri(c-1)) + (((c-(c-1))/((c+1)-(c-1)))*Ri(c+1));
c=c+2;
end
%realocate values to different matrix for influenza
S_lin_int=Si;
I_lin_int=Ii;
R_lin_int=Ri;
%calculate error
%S error
V_int_S= Si(1, 3:2:end); %interpolated values
V_model_S= S_influenza(1,3:2:end);
num_lin_S= (sum(V_int_S-V_model_S))^2; %numerator of error formula
E_lin_S= sqrt(num_lin_S/49); %compute error for S values for linear interpolation
%I error
V_int_I= Ii(1, 3:2:end); %interpolated values
V_model_I= I_influenza(1,3:2:end);
num_lin_I= (sum(V_int_I-V_model_I))^2; %numerator of error formula
E_lin_I= sqrt(num_lin_I/49); %compute error for S values for linear interpolation
%R error
V_int_R= Ri(1, 3:2:end); %interpolated values
V_model_R= R_influenza(1,3:2:end);
num_lin_R= (sum(V_int_R-V_model_R))^2; %numerator of error formula
E_lin_R= sqrt(num_lin_R/49); %compute error for S values for linear interpolation

%quadratic interpolation
Si_q=S;
Ii_q= I;
Ri_q= R;
c=3; %counter for interpolation
while c<= 97
Si_q(c)=((((c-(c+1))*(c-(c+3)))/(((c-1)-(c+1))*((c-1)-(c+3))))*Si_q(c-1)) + ((((c-(c-1))*(c-(c+3)))/(((c+1)-(c-1))*((c+1)-(c+3)))*Si_q(c+1))) + (((c-(c-1))*(c-(c+1)))/(((c+3)-(c-1))*((c+3)-(c+1))))*Si_q(c+3); %quadratic interpolation lagrange form
Ii_q(c)= ((((c-(c+1))*(c-(c+3)))/(((c-1)-(c+1))*((c-1)-(c+3))))*Ii_q(c-1)) + ((((c-(c-1))*(c-(c+3)))/(((c+1)-(c-1))*((c+1)-(c+3)))*Ii_q(c+1))) + (((c-(c-1))*(c-(c+1)))/(((c+3)-(c-1))*((c+3)-(c+1))))*Ii_q(c+3);
Ri_q(c)= ((((c-(c+1))*(c-(c+3)))/(((c-1)-(c+1))*((c-1)-(c+3))))*Ri_q(c-1)) + ((((c-(c-1))*(c-(c+3)))/(((c+1)-(c-1))*((c+1)-(c+3)))*Ri_q(c+1))) + (((c-(c-1))*(c-(c+1)))/(((c+3)-(c-1))*((c+3)-(c+1))))*Ri_q(c+3);
c=c+2;
end
%interpolate for the 99th value
Si_q(99)= (((99-95)*(99-97))/((93-95)*(93-97)))*Si_q(93) + (((99-93)*(99-97))/((95-93)*(95-97)))*Si_q(95) + (((99-93)*(99-95))/((97-93)*(97-95)))*Si_q(97); 
Ii_q(99)= (((99-95)*(99-97))/((93-95)*(93-97)))*Ii_q(93) + (((99-93)*(99-97))/((95-93)*(95-97)))*Ii_q(95) + (((99-93)*(99-95))/((97-93)*(97-95)))*Ii_q(97);
Ri_q(99)= (((99-95)*(99-97))/((93-95)*(93-97)))*Ri_q(93) + (((99-93)*(99-97))/((95-93)*(95-97)))*Ri_q(95) + (((99-93)*(99-95))/((97-93)*(97-95)))*Ri_q(97);
%realocate values to different matrix for influenza
S_quad_int=Si_q;
I_quad_int=Ii_q;
R_quad_int=Ri_q;

%calculate error
%S error
V_int_S= Si_q(1, 3:2:end); %interpolated values
V_model_S= S_influenza(1,3:2:end);
num_quad_S= (sum(V_int_S-V_model_S))^2; %numerator of error formula
E_quad_S= sqrt(num_quad_S/49); %compute error for S values for linear interpolation
%I error
V_int_I= Ii_q(1, 3:2:end); %interpolated values
V_model_I= I_influenza(1,3:2:end);
num_quad_I= (sum(V_int_I-V_model_I))^2; %numerator of error formula
E_quad_I= sqrt(num_quad_I/49); %compute error for S values for linear interpolation
% R error
V_int_R= Ri_q(1, 3:2:end); %interpolated values
V_model_R= R_influenza(1,3:2:end);
num_quad_R= (sum(V_int_R-V_model_R))^2; %numerator of error formula
E_quad_R= sqrt(num_quad_R/49); %compute error for S values for linear interpolation

%create a table to show the error values
Scol= [E_lin_S, E_quad_S]; %columb one
Icol= [E_lin_I, E_quad_I]; %columb two
Rcol= [E_lin_R, E_quad_R]; %columb three
Table= table(Scol', Icol', Rcol', 'VariableNames',{'Susceptible population','Infected population','Recovered population'}, 'RowNames',{'Linear Interpolation error','Quadratic Interpolation error'});
disp(Table)

