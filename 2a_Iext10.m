%made by jnullj - 2017
%huxley hodgkin nerve cell model
%membrane voltage for diff currents

%timing and stepping
dt = 0.01;		%step size
wtime = 100;		%total time to simulate
steps = wtime/dt;	%number of steps in time
%build empty arrays for data
V = zeros(1,steps);
n = zeros(1,steps);
m = zeros(1,steps);
h = zeros(1,steps);
%intial values
V(1) = 0;
n(1) = 0.3177;
m(1) = 0.0529;
h(1) = 0.5961;
c = 1;
g_K = 36;
g_Na = 120;
g_l = 0.3;
V_K = 12;
V_Na = -115;
V_l = -10.6;
%input 1
%I = zeros(1,steps); case 1
%I = ones(1,steps)*-5; case 2
I = ones(1,steps)*-10; %case 3
%functions
alpha_n = @(V) 0.01*(V+10)/(exp((V+10)/10)-1);
beta_n = @(V) 0.125*exp(V/80);
alpha_m = @(V) 0.1*(V+25)/(exp((V+25)/10)-1);
beta_m = @(V) 4*exp(V/18);
alpha_h = @(V) 0.07*exp(V/20);
beta_h = @(V) 1/(exp((V+30)/10)+1);
% v’=f(t,v)
f_v = (-g_K*(n(1)^4)*(V(1)-V_K)-g_Na*(m(1)^3)*h(1)*(V(1)-V_Na)-g_l*(V(1)-V_l)+I(1))/c;
f_n = alpha_n(V(1))*(1-n(1))-beta_n(V(1))*n(1);
f_m = alpha_m(V(1))*(1-m(1))-beta_m(V(1))*m(1);
f_h = alpha_h(V(1))*(1-h(1))-beta_h(V(1))*h(1);

%using euler method to find all values
for t=2:steps
	V(t) = V(t-1) + dt * f_v;
	n(t) = n(t-1) + dt * f_n;
	m(t) = m(t-1) + dt * f_m;
	h(t) = h(t-1) + dt * f_h;
	f_v = (-g_K*(n(t)^4)*(V(t)-V_K)-g_Na*(m(t)^3)*h(t)*(V(t)-V_Na)-g_l*(V(t)-V_l)+I(t))/c;
	f_n = alpha_n(V(t))*(1-n(t))-beta_n(V(t))*n(t);
	f_m = alpha_m(V(t))*(1-m(t))-beta_m(V(t))*m(t);
	f_h = alpha_h(V(t))*(1-h(t))-beta_h(V(t))*h(t);
end

V1 = -V-65;

%show plots
tval = 1:steps;

plot(tval, V1)
title('Iext=10')
xlabel('Time [ms]')
ylabel('Action potential [mV]')
ax = gca;
ax.XTick = 1:wtime:steps;
ax.XTickLabel = 0:wtime;
