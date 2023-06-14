% AEROSP-535
% Author: Jacob Dewey
% Final Project
%--------------------------------------------------------------------------
% This code simulates the launch of our designed vehicle from the pad into
% the desired orbit of 400km. It does so through a numerical simulation
% using a number of functions found at the bottom of the code. 
%--------------------------------------------------------------------------



clc, clear all, close all


%% Rocket Constants:
% Other Constants and Variables:
re = 6378000; % Radius of Earth [m]
g0 = 9.8; % Gravity of earth [m/s^2]
mu = 398600; %
ifpa = pi/2; % Initial Flight Angle [rad]
dt = 0.01; % Timestep [s]
targetAlt = 400e3; % Target Altitude [m]
altMargin = 500; % Margin for Altitude [m]
t_turn = 25; % Time of Gravity Turn [s]
init_dfpa = deg2rad(0.5); % Initial Gravity Turn Maneuver [rad]
orbSpeed = 7600; % Velocity Needed to Achieve Orbit [m/s]
speedMargin = 30; % Allowable Margin for Velocity [m/s]
fpaMargin = 15; % Allowable Margin for Flight Path Angle [deg]
check = 0;
check1 = 0;
check2 = 0;

%% Simulation:
% pitchAngleVec = deg2rad(linspace(20,35,1000));
pitchAngleVec = deg2rad(34);
pitchAngleVec2 = deg2rad(-8);
pitchj = 1;

for pitchi = 1:length(pitchAngleVec)
    pitchAngle = pitchAngleVec(pitchi);
    pitchAngle2 = pitchAngleVec2(pitchj);
    % Stage 1:
    s1.mi = 2873102.09; % Intial Stage Mass [kg]
    s1.tf = 79.1; % Burn time [s]
    s1.thrust = 34565000; % Stage Thrust [N]
    s1.mdot = 12588; % Mass Flow Rate of Engines [kg/s]
    s1.stage = 1;

    % Stage 2:
    s2.mi = 1825995.31; % Intial Stage Mass [kg]
    s2.tf = 396.5; % Burn time [s]
    s2.thrust = 18203400; % Stage Thrust [N]
    s2.mdot = 3650.08; % Mass Flow Rate of Engines [kg/s]
    s2.stage = 2;

    % Stage 3:
    s3.mi = 270157.33; % Intial Stage Mass [kg]
    s3.tf = 233.3; % Burn time [s]
    s3.thrust = 2022600; % Stage Thrust [N]
    s3.mdot = 521.44; % Mass Flow Rate of Engine [kg/s]
    s3.stage = 3;
    massEmpty = s3.mi - s3.mdot*s3.tf;

    % Intial Timestep:
    m = s1.mi - s1.mdot*dt;
    a = s1.thrust/m - g0*sin(ifpa);
    v = dt*a;
    alt = v*dt*sin(ifpa);
    xdis = v*dt*cos(ifpa);
    fpa = ifpa;
    theta = fpa;


    i = 2;
    for t = dt:dt:s1.tf+s2.tf+s3.tf
        if t < s1.tf-.00001 % First Stage Burn
            [m(i), a(i), v(i), alt(i), xdis(i), fpa(i), fpadot(i)] = sim(m(i-1), v(i-1), alt(i-1), xdis(i-1), fpa(i-1), fpa(i-1), pitchAngle2, dt, s1, t, t_turn, init_dfpa);
            if check1 == 0 && alt(i-1) > 40e3
                check1 = 1;
                pitch1 = i;
                pitch1t = t;
            elseif check2 == 0 && alt(i) > 200e3
                check2 = 1;
                pitch2 = i;
                pitch2t = t;
            end

            i = i+1;


        elseif abs(t-s1.tf) < 0.00001  % Staging Step from Stage 1 to Stage 2
            fprintf('Staging from 1 to 2: t + %1.2f seconds\n',t)
            meco = i-1;
            m(i) = s2.mi - s2.mdot*dt;
            a(i) = s2.thrust/m(i)*cos(fpa(i-1)-fpa(i-1)) - g0*sin(fpa(i-1));
            v(i) = v(i-1) + dt*a(i);
            alt(i) = alt(i-1) + v(i)*dt*sin(fpa(i-1));
            xdis(i) = xdis(i-1) + v(i)*dt*cos(fpa(i-1));
            r = re + alt(i);
            fpadot(i) = -(g0/v(i) - v(i)/r)*cos(fpa(i-1));
            fpa(i) = fpa(i-1) + dt*fpadot(i);
            gravityLoss(i) = g0*sin(fpa(i))*dt;

            if check1 == 0 && alt(i-1) > 40e3
                check1 = 1;
                pitch1 = i;
                pitch1t = t;
            elseif check2 == 0 && alt(i) > 200e3
                check2 = 1;
                pitch2 = i;
                pitch2t = t;
            end

            i = i+1;

        elseif t < s1.tf+s2.tf-.00001 % Second Stage Burn
            [m(i), a(i), v(i), alt(i), xdis(i), fpa(i), fpadot(i)] = sim(m(i-1), v(i-1), alt(i-1), xdis(i-1), fpa(i-1), pitchAngle, pitchAngle2, dt, s2, t, t_turn, init_dfpa);
            gravityLoss(i) = g0*sin(fpa(i))*dt;

            if check1 == 0 && alt(i-1) > 40e3
                check1 = 1;
                pitch1 = i;
                pitch1t = t;
            elseif check2 == 0 && alt(i) > 200e3
                check2 = 1;
                pitch2 = i;
                pitch2t = t;
            end

            i = i+1;

        elseif abs(t-s1.tf-s2.tf) < 0.00001 % Staging Step from Stage 2 to Stage 3
            seco = i-1;
            fprintf('Staging from 2 to 3: t + %1.2f seconds\n',t)
            if alt(i-1) >= 200e3
                pitchAngles3 = pitchAngle2;
            else
                pitchAngles3 = pitchAngle;
            end
            m(i) = s3.mi - s3.mdot*dt;
            a(i) = s3.thrust/m(i)*cos(pitchAngle-fpa(i-1)) - g0*sind(fpa(i-1));
            v(i) = v(i-1) + dt*a(i);
            alt(i) = alt(i-1) + v(i)*dt*sin(fpa(i-1));
            xdis(i) = xdis(i-1) + v(i)*dt*cos(fpa(i-1));
            r = re + alt(i);
            fpadot(i) = -(g0/v(i) - v(i)/r)*cos(fpa(i-1)) + s3.thrust/(v(i)*m(i))*sin(pitchAngles3-fpa(i-1));
            fpa(i) = fpa(i-1) + dt*fpadot(i);
            gravityLoss(i) = g0*sin(fpa(i))*dt;

            if check1 == 0 && alt(i-1) > 40e3
                check1 = 1;
                pitch1 = i;
                pitch1t = t;
            elseif check2 == 0 && alt(i) > 200e3
                check2 = 1;
                pitch2 = i;
                pitch2t = t;
            end

            i = i+1;

        elseif t < s1.tf+s2.tf+s3.tf-.00001 % Third Stage Burn
            [m(i), a(i), v(i), alt(i), xdis(i), fpa(i), fpadot(i)] = sim(m(i-1), v(i-1), alt(i-1), xdis(i-1), fpa(i-1), pitchAngle, pitchAngle2, dt,  s3, t, t_turn, init_dfpa);
            gravityLoss(i) = g0*sin(fpa(i))*dt;

            if check1 == 0 && alt(i-1) > 40e3
                check1 = 1;
                pitch1 = i;
                pitch1t = t;
            elseif check2 == 0 && alt(i) > 200e3
                check2 = 1;
                pitch2 = i;
                pitch2t = t;
            end

            i = i+1;
        end

        if abs(v(i-1)-orbSpeed) < speedMargin && abs(alt(i-1)-targetAlt) < altMargin && abs(rad2deg(fpa(i-1))) < fpaMargin
            teco = i-1;
            fprintf('Orbital Insertion Complete :) at t+%1.2f seconds\n\n',t)
            fprintf('Final Altitude: %1.3f km\n',alt(i-1)/1000)
            fprintf('Final Velocity: %1.3f km/s\n',v(i-1)/1000)
            fprintf('Final Mass: %1.2f kg \n',m(i-1))
            fprintf('Turn Time t+ %1.2f s\n',t_turn)
            fprintf('Turn Angle = %1.2f deg\n',rad2deg(init_dfpa))
            fprintf('Pitch Angle 1 = %1.2f deg at t+%1.2f seconds\n',rad2deg(pitchAngle),pitch1t)
            fprintf('Pitch Angle 2 = %1.2f deg at t+%1.2f seconds\n',rad2deg(pitchAngle2),pitch2t)
            check = 1;
            break
        elseif alt(i-1) < 0
            fprintf('Rocket Crashed :( at t+%1.2f seconds\n\n',t)
            break
        elseif abs(v(i-1)-orbSpeed) < speedMargin*0.1
            fprintf('Orbital Speed reached at t+%1.2f seconds\n',t)
            fprintf('Altitude = %1.2f km\n\n',alt(i-1)/1000)
            break
        elseif t > s1.tf+s2.tf+s3.tf-.00001
            fprintf('No Orbit, Ran Out of Fuel :( at t+%1.2f seconds\n',t')
            fprintf('Alt = %1.2f km, V = %1.2f km/s\n\n', alt(i-1)/1000, v(i-1)/1000)
            break
        end
    end
    
    if check == 1
        break
    end
end


% figure;
% plot(x,rad2deg(finalFPA))
% ylabel('Final FPA [deg]')
%
% figure;
% plot(x,finalVel/1000)
% ylabel('Final Vel [km/s]')
%
% figure;
% plot(x,finalAlt/1000)
% ylabel('Final Alt[km]')


%% Plots:
tl = length(alt);
tv = linspace(0,tl*dt,tl);
mecotxt = '\leftarrow MECO';
secotxt = '\leftarrow SECO';
tecotxt = 'Third Engine Cutoff, Orbit Achieved \rightarrow';
pitch1txt = '\leftarrow Pitch Manuever 1';
pitch2txt = '\leftarrow Pitch Manuever 2';


% Altitude vs Time Plot:
figure;
plot(tv,alt/1000,'LineWidth',2)
text(tv(meco),alt(meco)/1000,mecotxt)
text(tv(seco),alt(seco)/1000,secotxt)
text(tv(i-1),alt(i-1)/1000,tecotxt,'HorizontalAlignment','right')
text(tv(pitch1),alt(pitch1)/1000,pitch1txt)
text(tv(pitch2),alt(pitch2)/1000,pitch2txt)
title('Altitude vs Time')
xlabel('Time [sec]')
ylabel('Altitude [km]')
print('AltvsTime','-dpng','-r300')

% Mass vs Time Plot:
figure;
plot(tv,m,'LineWidth',2)
title('Mass vs Time')
xlabel('Time [sec]')
ylabel('Mass [kg]')
print('MassvsTime','-dpng','-r300')

% Altitude vs Downrange Distance Plot:
figure;
hold on
plot(xdis/1000,alt/1000,'LineWidth',2)
text(xdis(meco)/1000,alt(meco)/1000,mecotxt)
text(xdis(seco)/1000,alt(seco)/1000,secotxt)
text(xdis(i-1)/1000,alt(i-1)/1000,tecotxt,'HorizontalAlignment','right')
text(xdis(pitch1)/1000,alt(pitch1)/1000,pitch1txt)
text(xdis(pitch2)/1000,alt(pitch2)/1000,pitch2txt)
title('Altitude vs Downrange Distance')
xlabel('Downrange Distance [km]')
ylabel('Altitude [km]')
print('AltvsXdis','-dpng','-r300')

% Acceleration vs Time Plot:
figure;
plot(tv,a,'LineWidth',2)
title('Acceleration vs Time')
xlabel('Time [sec]')
ylabel('Acceleration [m/s^2]')
print('AccvsTime','-dpng','-r300')

% G-Loading vs Altitude Plot:

figure;
plot(alt/1000,a/g0,'LineWidth',2)
text(alt(meco)/1000,a(meco)/g0,mecotxt)
text(alt(seco)/1000,a(seco)/g0,secotxt)
text(alt(i-1)/1000,a(i-1)/g0,tecotxt,'HorizontalAlignment','right')
text(alt(pitch1)/1000,a(pitch1)/g0,pitch1txt)
text(alt(pitch2)/1000,a(pitch2)/g0,pitch2txt)
title('G-Loading vs Altitude')
xlabel('Altitude [km]')
ylabel('G-Loading')
print('GvsTime','-dpng','-r300')

% Velocity vs Time Plot:
figure;
plot(tv,v,'LineWidth',2)
title('Velocity vs Time')
xlabel('Time [sec]')
ylabel('Velocity [m/s]')
print('VelvsTime','-dpng','-r300')

% FPA vs Time Plot:
figure;
plot(tv,rad2deg(fpa),'LineWidth',2)
title('Flight Path Angle vs Time')
xlabel('Time [sec]')
ylabel('Flight Path Angle [deg]')
print('FPAvsTime','-dpng','-r300')


%
% figure;
% plot(tv(20.51/dt:tl),rad2deg(fpadot(20.51/dt:tl)),'LineWidth',2)
% title('dfpa vs Time')
% xlabel('Time [sec]')
% ylabel('dfpa [deg/s]')

%% Functions:
function [m, a, v, alt, xdis, fpa, fpadot] = sim(mIn, vIn, altIn, xdisIn, fpaIn, thetaIn, thetaIn2, dt, stage, t, kickTime, kickInit)


if t < kickTime
    re = 6378000;
    g0 = 9.8; % Gravity of earth
    m = mIn - stage.mdot*dt;
    a = (stage.thrust/m) - g0*sin(fpaIn);
    v = vIn + a*dt;
    alt = altIn + v*sin(fpaIn)*dt;
    xdis = xdisIn + v*cos(fpaIn)*dt;
    r = re+alt;
    fpadot = 0;
    fpa = fpaIn;
elseif abs(t-kickTime) < 0.000001
    re = 6378000;
    g0 = 9.8; % Gravity of earth
    m = mIn - stage.mdot*dt;
    a = (stage.thrust/m) - g0*sin(fpaIn);
    v = vIn + a*dt;
    alt = altIn + v*sin(fpaIn)*dt;
    xdis = xdisIn + v*cos(fpaIn)*dt;
    r = re+alt;
    fpadot = -kickInit;
    fpa = fpaIn - kickInit;
elseif altIn < 40e3
    re = 6378000;
    g0 = 9.8; % Gravity of earth
    m = mIn - stage.mdot*dt;
    a = (stage.thrust/m) - g0*sin(fpaIn);
    v = vIn + a*dt;
    alt = altIn + v*sin(fpaIn)*dt;
    xdis = xdisIn + v*cos(fpaIn)*dt;
    r = re+alt;
    fpadot = - (g0/v - v/(r))*cos(fpaIn);
    fpa = fpaIn + dt*fpadot;
elseif altIn >= 40e3 && altIn < 200e3
    re = 6378000;
    g0 = 9.8; % Gravity of earth
    m = mIn - stage.mdot*dt;
    a = (stage.thrust/m*cos(thetaIn-fpaIn)) - g0*sin(fpaIn);
    v = vIn + a*dt;
    alt = altIn + v*sin(fpaIn)*dt;
    xdis = xdisIn + v*cos(fpaIn)*dt;
    r = re+alt;
    fpadot = - (g0/v - v/r)*cos(fpaIn) + stage.thrust/(v*m)*sin(thetaIn-fpaIn);
    fpa = fpaIn + dt*fpadot;
elseif altIn >= 200e3
    re = 6378000;
    g0 = 9.8; % Gravity of earth
    m = mIn - stage.mdot*dt;
    a = (stage.thrust/m*cos(thetaIn-fpaIn)) - g0*sin(fpaIn);
    v = vIn + a*dt;
    alt = altIn + v*sin(fpaIn)*dt;
    xdis = xdisIn + v*cos(fpaIn)*dt;
    r = re+alt;
    fpadot = - (g0/v - v/r)*cos(fpaIn) + stage.thrust/(v*m)*sin(thetaIn2-fpaIn);
    fpa = fpaIn + dt*fpadot;
end

end

