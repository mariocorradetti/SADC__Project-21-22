clc 
close all
clearvars
addpath('functions\')
addpath('functions\time\')

%% Orbit Data

% Set useful data
data.mu = astroConstants(13);  % Earth’s gravitational parameter [km^3/s^2]
data.we = deg2rad(15.04/3600); % Earth’s rotation velocity [km/s]
data.J2 = astroConstants(9);   % Gravitatonal field constant of the Earth  
data.Re = astroConstants(23);  % Earth's mean radius [km]
data.th_g0 = 0;                % Longitude of Greenwich meridian at initial time [rad]

% Given orbit 
a = 8150;   % Semi-major axis [km]
e = 0.0054;     % Eccentricity [-]
i =deg2rad(7.9875); % Inclination [rad]
OM =deg2rad(64.26);  % Our choice [rad]
om =deg2rad(348.04); % Our choice [rad]
th =deg2rad(0);  % Our choice [rad]
orbitpar0 = [a e i OM om th];
T = 2*pi*sqrt(a^3/data.mu); % Period of one orbit [s]
n = sqrt(data.mu/a^3); % Mean angular velocity

% Compute initial state vector
kep = [a,e,i,OM,om,th];
[rr0,vv0] = kep2car(kep, data.mu);
n0 = sqrt(data.mu / norm(rr0)^3);

% Apparent Sun orbit
T_sun = 1*365*24*60*60;
R_sun = astroConstants(2);
epsilon =  deg2rad(23.45);
n_sun = 2*pi / T_sun;

%% Enviroment 
% Sun radiation
Fe = astroConstants(31); %[W/m^2]
c = astroConstants(5);%[km/s]
c = c*1e3;
Fe_tot = Fe+400+(117+89)/2;
P_srp = Fe_tot/c;

%Magnetic
H0 = sqrt((-29615e-9)^2+(-1728e-9)^2+(5186e-9)^2);
R3H0 = data.Re^3*H0;

%% SENSORS

% MAGNETOMETER (NanoSense M315)
magnetometer.rate = 10; % [Hz] up to 140 
magnetometer.samplingtime = 1/magnetometer.rate;
magnetometer.maxval = 800e-6; %[T]
magnetometer.minval = -800e-6; %[T]

a_x = [1                    deg2rad(rand(1))            deg2rad(rand(1))];
a_y = [deg2rad(rand(1))              1                  deg2rad(rand(1))];
a_z = [deg2rad(rand(1))     deg2rad(rand(1))                    1       ];
magnetometer.A_eps = [a_x;a_y;a_z];
magnetometer.noise = 15e-9*[1 1 1]; 

% SUN SENSOR (FSS100 Nano Fine Sun Sensor)
sunsensor.FOV = 60 ; %[deg] (+- 60°)
sunsensor.accuracy = 0.1; %[deg] at 1 sigma (+-0.1°) STD with mean 0
sunsensor.rate = 8; %[Hz]
sunsensor.samplingtime = 1/sunsensor.rate; 

% EARTH HORIZON SENSOR (CubeSense N)
horizonsensor.FOV = 90; %[deg] (+-90°)
horizonsensor.rate = 2; %[Hz]
horizonsensor.samplingtime = 1 / horizonsensor.rate;
horizonsensor.accuracy = 0.2; % < 0.2° 3-sigma

weights = [0.15;0.5;0.35];

% GYROSCOPE 
gyro.rate = 20 ;% 250 [samples/sec] 20 to speed up the simulation 
gyro.samplingtime = 1/gyro.rate;
gyro.bias = 0.3; %[°/h]
gyro.ARW = 0.15; %[°/sqrt(h)]


%% Actuators
% Reaction wheel (BTC Micro Reaction wheel)
wheel.A = eye(3);
wheel.maxtorque = 0.6e-3; %[Nm]
wheel.momentum = 18e-3; %[Nms]
%Inertia wheel
wheel.maxrotationrate =6000; %[rpm] 
wheel.w = [0 0 6000*2*pi/60];
wheel.mass = 150e-3; %[kg]
wheel.R = (43e-3)/2; %[m]
wheel.R_in = (43e-3)/2; %[m]
wheel.h = 18e-3; %[m]
wheel.V_in=pi*wheel.R_in^2*wheel.h; %[m^3]
wheel.mass_in = 150e-3; %[kg]
wheel.rho=wheel.mass/wheel.V_in;
%off-design
wheel.R_in = (43e-3)/2; %[m]
wheel.mass_in = wheel.rho*(pi*wheel.R_in^2*wheel.h);
%
wheel.Ix = 1/12*wheel.mass_in*wheel.h^2 + 1/4*wheel.mass_in*wheel.R_in^2;
wheel.Iy = 1/12*wheel.mass_in*wheel.h^2 + 1/4*wheel.mass_in*wheel.R_in^2;
wheel.Iz = 1/2*wheel.mass_in*wheel.R_in^2;
wheel.I = diag([wheel.Ix wheel.Iy wheel.Iz]);

%% S/C Data
%{
                                      ^
                                      |
                                      |   x
                                      |
                             +--------|---------+
                            /         |        /|
                           /                  / |
                          /                  /  |
                         +------------------+   |  
                         |                  |   |
                         |                  |   |
                         |                  |   |
                         |                  |   |
                         |                  |   |h
               y <------ |                  |   |
                         |                  |   |
                         |                  |   |
                         |                  |   |
                         |                  |   |
                         |                  |  / 
                         |        w         | /  d
                         +------------------+
                                  /
                                 /  z
                                \/
%}
% Initial conditions
wx0 = 0.1;
wy0 = 0.1;
wz0 = 0.1;
w0 = [wx0 wy0 wz0]'; % Initial angular velocity 
Anb0 = eye(3);
Abn0 = Anb0'; % Inertial to body

% Assuming the s/c as an homogeneous solid 
sc.m_sc =3.5; % [kg]
sc.w = 0.1; % [m]
sc.d = 0.1; % [m]
sc.h = 0.3405; %[m]
sc.msp = 0.270; %[kg]
sc.th = 1.6e-3; %[m]

sc.m = sc.m_sc + wheel.mass*3 + wheel.mass_in;
sc.rcg = ([0 0 0]*sc.m + [sc.h/2-sc.th/2 sc.w/2+sc.h/2 0]*sc.msp/2 + [sc.h/2-sc.th/2 -(sc.w/2+sc.h/2) 0]*sc.msp/2)/(sc.m+sc.msp);

sc.Iz = 1/12*sc.m*(sc.h^2+sc.w^2);
sc.Iy = 1/12*sc.m*(sc.d^2+sc.h^2);
sc.Ix = 1/12*sc.m*(sc.w^2+sc.d^2);

sc.Ix_sp = 1/12*sc.msp/2*(sc.d^2+sc.h^2);
sc.Iy_sp = 1/12*sc.msp/2*(sc.d^2+sc.th^2);
sc.Iz_sp = 1/12*sc.msp/2*(sc.h^2+sc.th^2);
sc.rcg_sp1 = [sc.h/2 sc.w/2+sc.h/2 0]';
sc.rcg_sp2 = [sc.h/2 -(sc.w/2+sc.h/2) 0]';

sc.I_sp1 = diag([sc.Ix_sp sc.Iy_sp sc.Iz_sp])+sc.msp/2*( norm(sc.rcg_sp1)*eye(3)*kron(sc.rcg_sp1,sc.rcg_sp1'));
sc.I_sp2 = diag([sc.Ix_sp sc.Iy_sp sc.Iz_sp])+sc.msp/2*( norm(sc.rcg_sp2)*eye(3)*kron(sc.rcg_sp2,sc.rcg_sp2'));
sc.I = diag([sc.Ix,sc.Iy,sc.Iz])+sc.I_sp1+sc.I_sp2;
sc.Iinv = inv(sc.I);

% S/C main body
A_x = sc.w*sc.h;
A_y = sc.d*sc.h;
A_z = sc.d*sc.w;
A_sc = [A_x A_y A_z A_x A_y A_z];
% s/c main body normal surfaces
n1 = [1;0;0];
n2 = [0;1;0];
n3 = [0;0;1];
n4 = -n1;
n5 = -n2;
n6 = -n3;
% s/c main body surfaces parameter
rho_s_sc = 0.1;
rho_d_sc = 0.5;

% Solar panels
A_sp = (sc.w*sc.h).*ones(1,4);

n7 = [1;0;0];
n8 = -n7;
n9 = [1;0;0];
n10 = -n9;

rho_s_sp = 0.1;
rho_d_sp = 0.1;

% Distance from the CG
N = [n1 n2 n3 n4 n5 n6 n7 n8 n9 n10];

% center face to CG
r_sc_x = sc.d/2;
r_sc_y = sc.w/2;
r_sc_z = sc.h/2;
ri_sc = [r_sc_x , r_sc_y , r_sc_z , r_sc_x , r_sc_y , r_sc_z];
ri_sc = ri_sc.*N(:,1:6);
% center solar panel to CG
r_sp1 = [sc.h/2 0 sc.w/2+sc.h/2]';
r_sp2 = [sc.h/2 0 -sc.w/2-sc.h/2]';

ri = [ri_sc,r_sp1,r_sp1,r_sp2,r_sp2];

rho_s = cat(2,ones(1,6)*rho_s_sc,ones(1,4)*rho_s_sp);
rho_s(2)=rho_s_sp; % solar panels also on the s/c in y and -y direction
rho_s(5)=rho_s_sp;

rho_d = cat(2,ones(1,6)*rho_d_sc,ones(1,4)*rho_d_sp);

A = [A_sc,A_sp];% area vector

m = [0.01 0.05 0.01]'; %[A m^2]


%% From Simulink

open("Project_SADC.slx"); %for Matlab 2021b
% open("Project_SADC_2018.slx"); %for Matlab 2018b
out = sim('Project_SADC');

time = out.tout;
Abn = out.A_bn.Data;
Anl = out.A_nl.Data;
r_orbit = out.r_orbit.Data;
r_Sun = out.r_Sun_N.Data;
w = squeeze(out.w.Data);% Angular velocity
angerr = out.ang_err.Data;      % Pointing error 
M_SRP = squeeze(out.M_srp.Data);         % SRP Torque
M_GG = squeeze(out.M_GG.Data);           % Gravity Gradient Torque
M_mag = squeeze(out.M_magnetic.Data);    % Magnetic Torque
M_RW = out.torque_RW.Data;      % Reaction Wheel Torque

 
%% Plot

%Plot pointing error
figure(1)
plot(time(16191:end),angerr(16191:end),'LineWidth',1.5)
grid on
xlabel('Time [s]')
ylabel('Pointing Error [°]')
xlim([time(16191) time(end)])
ylim([-5 100])

%Angular Velocity
figure(2)
tiledlayout(2,3) 
ax1 = nexttile([1 3]);
plot(time,w(1,:),time,w(2,:),time,w(3,:),'LineWidth',1.5)
xlim([0 time(end)])
grid on
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')
legend('\omega_x','\omega_y','\omega_z')

ax2 = nexttile;
plot(time(1:16190),w(1,1:16190),time(1:16190),w(2,1:16190),time(1:16190),w(3,1:16190),'LineWidth',1.5)
xlim([time(1) time(16190)])
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')
title(ax2,'De-tumbling')
grid on

ax3 = nexttile;
plot(time(16191:80242),w(1,16191:80242),time(16191:80242),w(2,16191:80242),time(16191:80242),w(3,16191:80242),'LineWidth',1.5)
xlim([time(16191) time(80242)])
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')
title(ax3,'Slew maneuvre')
grid on

ax4 = nexttile;
plot(time(216250:end),w(1,216250:end),time(216250:end),w(2,216250:end),time(216250:end),w(3,216250:end),'LineWidth',1.5)
xlim([time(216250) time(end)])
xlabel('Time [s]')
ylabel('Angular velocity [rad/s]')
title(ax4,'Pointing')
grid on

% External Torque
figure(3)
plot(time,vecnorm(M_SRP),'LineWidth',1.5)
hold on
plot(time,vecnorm(M_GG),'LineWidth',1.5)
plot(time,vecnorm(M_mag),'LineWidth',1.5)
xlim([0 time(end)])
legend('M_{SRP}','M_{GG}','M_{mag}')
grid on
xlabel('Time [s]')
ylabel('External Torque [Nm]')

% Actuator torque
figure(4)
tiledlayout(2,3) 
ax1 = nexttile([1 3]);
plot(time,M_RW(:,1),time,M_RW(:,2),time,M_RW(:,3),'LineWidth',1.5)
xlim([0 time(end)])
grid on
xlabel('Time [s]')
ylabel('Torque [Nm]')
legend('T_x','T_y','T_z')

ax2 = nexttile;
plot(time(1:16190),M_RW(1:16190,1),time(1:16190),M_RW(1:16190,2),time(1:16190),M_RW(1:16190,3),'LineWidth',1.5)
xlim([time(1) time(16190)])
xlabel('Time [s]')
ylabel('Torque [Nm]')
title(ax2,'De-tumbling')
grid on

ax3 = nexttile;
plot(time(16160:80242),M_RW(16160:80242,1),time(16160:80242),M_RW(16160:80242,2),time(16160:80242),M_RW(16160:80242,3),'LineWidth',1.5)
xlim([time(16160) time(80242)])
xlabel('Time [s]')
ylabel('Torque [Nm]')
title(ax3,'Slew maneuvre')
grid on

ax4 = nexttile;
plot(time(216250:end),M_RW(216250:end,1),time(216250:end),M_RW(216250:end,2),time(216250:end),M_RW(216250:end,3),'LineWidth',1.5)
xlim([time(216250) time(end)])
xlabel('Time [s]')
ylabel('Torque[Nm]')
title(ax4,'Pointing')
grid on

%% Animated plot
% %uncomment
% figure('units','normalized','outerposition',[0 0 1 1])
% check=[1:200:length(time)];
% check=check(end);
% for i = 1:200:length(time)
% plot3(r_orbit(:,1)/1e3, r_orbit(:,2)/1e3, r_orbit(:,3)/1e3,'LineWidth', .5)
% % Inertial
% line([0 1],[0 0],[0 0],'Linewidth',2)
% grid on
% hold on
% line([0 0],[0 1],[0 0],'Linewidth',2)
% line([0 0],[0 0],[0 1],'Linewidth',2)
% axis(8*[-1.1 1.1 -1.1 1.1 -1.1 1.1])
% view(70.656353035073735,7.975234447350012)
% Anb = Abn(:,:,i)';
% % Body in inertial
% line(r_orbit(i,1)/1e3+[0 Anb(1,1)],r_orbit(i,2)/1e3+[0 Anb(2,1)],r_orbit(i,3)/1e3+[0 Anb(3,1)],'Color','k','Linewidth',2)
% line(r_orbit(i,1)/1e3+[0 Anb(1,2)],r_orbit(i,2)/1e3+[0 Anb(2,2)],r_orbit(i,3)/1e3+[0 Anb(3,2)],'Color','k','Linewidth',2)
% line(r_orbit(i,1)/1e3+[0 Anb(1,3)],r_orbit(i,2)/1e3+[0 Anb(2,3)],r_orbit(i,3)/1e3+[0 Anb(3,3)],'Color','k','Linewidth',2)
% % LVLH in inertial
% line(2*[0 Anl(1,1,i)],2*[0 Anl(2,1,i)],2*[0 Anl(3,1,i)],'Color','r','Linewidth',1)
% line(2*[0 Anl(1,2,i)],2*[0 Anl(2,2,i)],2*[0 Anl(3,2,i)],'Color','r','Linewidth',1)
% line(2*[0 Anl(1,3,i)],2*[0 Anl(2,3,i)],2*[0 Anl(3,3,i)],'Color','r','Linewidth',1)
% plot3([0 r_orbit(i,1)/1e3],[0 r_orbit(i,2)/1e3],[0 r_orbit(i,3)/1e3],'--')
% % s/c
% sat(Anb,r_orbit(i,:)/1e3,5)
% Plot_Earth(0,0,0)
% % sat(Anb,[0 0 0],10)
% axis([-10 10 -10 10 -10 10]);
% 
% pause(0.001)
% 
% if i<check
%     clf
% end
% 
% end
% 
