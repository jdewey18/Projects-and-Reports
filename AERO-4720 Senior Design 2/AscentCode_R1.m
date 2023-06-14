% AERO-4720
% Author: Jacob Dewe
% Calculation of speed, mach number, and altitude:
clc, clear all, close all 
format short g

%% Force Calculations:
% Constants
[T, P, rho] = atmo_model_100km(100); % Atmospheric model (to 100km)
Cd = 0.27; % Drag Coefficient: best guess
s = 72.63; % Planform Surface Area [m^2]
g0 = 9.8067; % Gravitational acceleration [m/s^2]
mu = 398600; %
gam = pi/2;  % Initial flight path angle (measured from the horizontal) [deg]
mass_0 = 3076259.09; % Zero Stage (SRBs) Initial Mass [kg]
mass_1 = 1944415.29; % First Stage Initial Mass [kg]
mass_2 = 325559.02;  % Second Stage Initial Mass [kg]
burnoutMass1 = 431673.29; % First Stage "Dry" Mass [kg]
burnoutMass2 = 98530.96;  % Second Stage "Dry" Mass [kg]
minThrottle1 = 0.65; % Minimum Throttle on First Stage Engines
minThrottle2 = 0.39; % Minimum Throttle on Second Stage Engines
fuelmass_0 = 1092000; % Zero Stage Fuel Mass [kg]
fuelmass_1 = 1512472; % First Stage Fuel Mass [kg]
fuelmass_2 = 227028;  % Second Stage Fuel Mass [kg]
T0 = 28.32e6; % Zero Stage Thrust [N]
T1 = 9.6e6;   % First Stage Thrust [N]
T2 = 480e3;   % Second Stage Thrust [N]
mdot_0 = 8000;    % Solid Booster Fuel Mass Flow Rate [kg/s]
mdot_1 = 3272.88; % First Stage Fuel Mass Flow Rate [kg/s]
mdot_2 = 157.33;  % Second Stage Fuel Mass Flow Rate [kg/s]
re = 6378e3; % Radius of Earth [m]
gamInit = 1; % Initial kickover angle in degrees
deg1 = 0.5;
deg2 = 0.7;
deg3 = 1.6;
dt = 0.01; % Timestep [s]
gTurnTime = 40/dt; % Number of iterations until gravity turn [seconds/seconds]

% Calculations
burnTime0 = fuelmass_0/(mdot_0); % Zero Stage Burn Time [s]
burnTime1 = fuelmass_1/(mdot_1); % First Stage Burn Time [s]
burnTime2 = fuelmass_2/(mdot_2); % Second Stage Burn Time [s]
totalBurn = burnTime1 + burnTime2; % Total Burn Time of Vehicle [s]


%% Ascent Calculations:
% Initial:
burnmass = mdot_0+mdot_1; % Initial Fuel Mass Consumption [kg/s]
fmass(1) = mass_0 - dt*burnmass; % Stage Fuel Mass Remaining [kg]
D(1) = 0; % Drag [N]
acc(1) = (T0 + T1)/fmass -D(1)/fmass(1) - g0*sin(gam); % Acceleration of Vehicle [m/s^2]
% racc(1) = (T0+T1)/fmass; % Relative Acceleration (used for code testing, not an actual variable) [m/s^2]
v(1) = dt*acc(1); % Velocity of Vehicle [m/s]
alt(1) = dt*v(1)*sin(gam); % Altitude of Vehicle [m]
xdis(1) = dt*v(1)*cos(gam); % Downrange Distance of Vehicle [m]
gam(1) = gam; % Flight Path Angle of Vehicle [rad]

% Loop
orbit = false; % Sets Orbit Variable [boolean]
t = dt; % Sets Initial time [s]
i = 2; % Iteration Number
while fmass(i-1) > burnoutMass2 % Checks that the stage has fuel remaining
    if t(i-1) < burnTime0 % First Stage burn with boosters
        g(i) = g0/(1+alt(i-1)/re)^2; % Gravitational Acceleration Felt by Vehicle [m/s^2]
        if alt/1000 < 100 % Checks to see if still in atmosphere (below 100km)
            D(i) = dragCalc(v(i-1), rho(fix(alt(i-1)*1e-3)+1), Cd, s);
        else 
            D(i) = 0;
        end
        [acc(i), fmass(i), racc(i)] = accFunction(fmass(i-1), (T0+T1), (mdot_1+mdot_0), D(i), 1, g(i), gam(i-1), dt, (t(i-1)+dt));
        [v(i), alt(i), xdis(i)] = trajectoryFunction(alt(i-1), xdis(i-1), v(i-1), acc(i), dt, gam(i-1));
        if i < gTurnTime
            gam(i) = gam(i-1);
        elseif i == gTurnTime
            gam(i) = gam(i-1) - gamInit*pi/180;
        else
            gam(i) = gamCalc(gam(i-1), g(i), v(i-1), alt(i), deg1, dt);
        end
        
    elseif  abs(t(i-1)+dt - burnTime0) < 0.001 % Staging the boosters     t(i-1) == burnTime0 
        g(i) = g0/(1+alt(i-1)/re)^2;
        if alt/1000 < 100
            D(i) = dragCalc(v(i-1), rho(fix(alt(i-1)*1e-3)+1), Cd, s);
        else 
            D(i) = 0;
        end
        [acc(i), fmass(i), racc(i)] = accFunction((mass_1 - i*dt*mdot_1), (T1), (mdot_1), D(i), 1, g(i), gam(i-1), dt, (t(i-1)+dt));
        [v(i), alt(i), xdis(i)] = trajectoryFunction(alt(i-1), xdis(i-1), v(i-1), acc(i), dt, gam(i-1));
        gam(i) = gamCalc(gam(i-1), g(i), v(i-1), alt(i), deg1, dt);
        
        % fprintf('Booster Separation at t+%3.2fs \n\n', t(i-1)+dt)
        
    elseif t(i-1) > burnTime0 && fmass(i-1) > (burnoutMass1 + mdot_1) % Rest of First Stage Burn   
        g(i) = g0/(1+alt(i-1)/re)^2;
        if alt/1000 < 100
            D(i) = dragCalc(v(i-1), rho(fix(alt(i-1)*1e-3)+1), Cd, s);
        else 
            D(i) = 0;
        end
        [acc(i), fmass(i), racc(i)] = accFunction(fmass(i-1), (T1), (mdot_1), D(i), minThrottle1, g(i), gam(i-1), dt, (t(i-1)+dt));
        [v(i), alt(i), xdis(i)] = trajectoryFunction(alt(i-1), xdis(i-1), v(i-1), acc(i), dt, gam(i-1));
        gam(i) = gamCalc(gam(i-1), g(i), v(i-1), alt(i), deg2, dt);
        
        
    elseif t(i-1) > burnTime0 && fmass(i-1) > burnoutMass1 && fmass(i-1) < (burnoutMass1 + mdot_1) % Staging between First and Second Stage
        g(i) = g0/(1+alt(i-1)/re)^2;
        if alt/1000 < 100
            D(i) = dragCalc(v(i-1), rho(fix(alt(i-1)*1e-3)+1), s, dt);
        else 
            D(i) = 0;
        end
        [acc(i), fmass(i), racc(i)] = accFunction(mass_2, (T2), (mdot_2), D(i), minThrottle2, g(i), gam(i-1), dt, (t(i-1)+dt));
        [v(i), alt(i), xdis(i)] = trajectoryFunction(alt(i-1), xdis(i-1), v(i-1), acc(i), dt, gam(i-1));
        gam(i) = gamCalc(gam(i-1), g(i), v(i-1), alt(i), deg3, dt);
        
        fprintf('First Stage Separation at t+%3.2fs \n\n', t(i-1)+dt)
        
        
    else % Second Stage Burn
        g(i) = g0/(1+alt(i-1)/re)^2;
        D(i) = 0;
        [acc(i), fmass(i), racc(i)] = accFunction(fmass(i-1), T2, (mdot_2), D(i), minThrottle2, g(i), gam(i-1), dt, (t(i-1)+dt));
        [v(i), alt(i), xdis(i)] = trajectoryFunction(alt(i-1), xdis(i-1), v(i-1), acc(i), dt, gam(i-1));
        gam(i) = gamCalc(gam(i-1), g(i), v(i-1), alt(i), deg3, dt);
        
        
    end
    
%     disp(acc(i))
%     disp(v(i))
    
    t(i) = t(i-1) + dt;
    
    dragLoss(i) = D(i)*dt/fmass(i);
    gravityLoss(i) = g(i)*sin(gam(i))*dt;
    
    vOrbit = sqrt(mu*1000/(re+alt(i)));
    if v(i)/1000 + 0.40859 - vOrbit > 0 && abs(gam(i)) < 0.001
        deltaV = 3050.91*log(fmass(i)/burnoutMass2);
        orbit = true;
        disp('Success! You''re in orbit!')
        fprintf('t+%6.2f seconds\n', t(i))
        fprintf('Fuel Left = %9.2f kg\n', (fmass(i) - burnoutMass2))
        fprintf('Delta-V Left = %4.2f m/s\n\n', deltaV)
        break
    end
    
    
    if alt(i) < 0
        disp('You crashed :(')
        break
    end
        
    
    i = i+1;
end

if orbit == true
    totalDrag = sum(dragLoss);
    totalGravity = sum(gravityLoss);
    targetR = (re+800000)/1000;
    vTarget = sqrt(mu/targetR);
    maneuverV = abs(v(length(v))+408.59 - vTarget*1000);
    
    fprintf('Drag Losses: %4.2f m/s\n', totalDrag)
    fprintf('Gravity Losses: %4.2f m/s\n', totalGravity)
    fprintf('Delta-V to Target: %3.2f m/s \n', maneuverV)
    fprintf('Final Delta-V: %3.2f m/s \n\n', deltaV - maneuverV)
    
end

%% Plots:
% Alt vs Time:
figure;
plot([1:1:length(alt)]*dt,alt/1000,'Linewidth',2)
title('Altitude vs Time')
xlabel('Time [s]')
ylabel('Altitude [km]')

% Acceleration vs Time
figure;
plot([1:1:length(acc)]*dt,acc,'Linewidth',2) 
title('Acceleration vs Time')
xlabel('Time [s]')
ylabel('Accelaration [m/s^2]')

% Velocity vs Time
figure;
plot([1:1:length(v)]*dt,v,'Linewidth',2) 
title('Velocity vs Time')
xlabel('Time [s]')
ylabel('Velocity [m/s]')

% % No-Loss Acceleration vs Time
% figure;
% plot([1:1:length(racc)]*dt,racc,'Linewidth',2) 
% title('Reference Acceleration vs Time')
% xlabel('Time [s]')
% ylabel('Acceleration [m/s^2]')

% Drag vs Time
figure;
plot([1:1:length(D)]*dt,D,'Linewidth',2)
title('Drag vs Time')
xlabel('Time [s]')
ylabel('Drag [N]')

% Mass vs Time
figure;
plot([1:1:length(fmass)]*dt,fmass,'Linewidth',2)
title('Mass vs Time')
xlabel('Time [s]')
ylabel('Mass [kg]')


% for i = 1:length(v)
%     gdr(i) = -(g(i)/v(i) - v(i)/(alt(i)+re));
% end
% 
% figure;
% plot([1:length(gdr)]*dt,gdr)


% FPA vs Time
figure;
plot([1:length(gam)]*dt,gam*180/pi)
title('Flight Path Angle vs Time')
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')

% FPA vs Altitude
figure;
plot(alt/1000,gam*180/pi)
title('Flight Path Angle vs Altitude')
xlabel('Alt [km]')
ylabel('Flight Path Angle [deg]')

% 
% figure;
% plot([1:length(g)]*dt,g)



%% Functions:
% Gravity Turn:
function gam = gamCalc(gamIn, g, v, alt, deg, dt)
real = 0;
if real == 1
    Re = 6378*1000; % Radius of the Earth [km]
    r = Re+alt;
    gDot = -(g/v - v/r)*cos(gamIn);
    gam = gamIn + gDot*dt;
elseif real == 0
    gDot = -deg*sin(gamIn)*pi/180;
    gam = gamIn + gDot*dt;
end
end

% Drag Calculator
function D = dragCalc(v, rho, Cd, s)
D = 0.5*Cd*s*rho*v^2;

end

% Acceleration Function (includes a throttling mechanism):
function [acc, fMass, rAcc] = accFunction(fMassIn, thrust, mDotIn, D, minThrottle, g, gam, dt, t)
g0 = 9.8067; % m/s^2
refAcc = (thrust - D)/(fMassIn-mDotIn) - g*sin(gam);
rAcc = thrust/(fMassIn-mDotIn);

% Ensures the acceleration is under 2 g's
if refAcc <= 3*g0  % If under 2 g's, no throttling occurs
    acc = refAcc;
    fMass = fMassIn - dt*mDotIn;
elseif refAcc > 3*g0  % If greater than 2 g's, it throttles the engines
    disp(refAcc)
    desThrust = (3*g0 + g*sin(gam))*fMassIn + D;
    throttle = desThrust/thrust;
    if throttle > minThrottle   % Ensures the throttle doesn't exceed limits
        mDot = dt*throttle*mDotIn;
        fMass = fMassIn - mDot;
        acc = (throttle*thrust - D)/fMass - g*sin(gam);
        fprintf('Throttled to %1.2f%% at t+%1.2fs \n\n', throttle*100, t)
    else  % If throttle ability is exceeded, ensures lowest g's possible
        mDot = dt*minThrottle*mDotIn;
        fMass = fMassIn - mDot;
        acc = (minThrottle*thrust - D)/fMass - g*sin(gam);
        fprintf('Throttled to %1.2f%% at t+%1.2fs \n\n', minThrottle*100, t)
    end
end
end

% Trajector Finder:
function [v, alt, xdis] = trajectoryFunction(altIn, xdisIn, vIn, accIn, dt, gam)
v = vIn + dt*accIn;             % Determines the velocity [m/s]
alt = altIn + dt*v*sin(gam);    % Determines the altitude [m]
xdis = xdisIn + dt*v*cos(gam);  % Determines downrange distance [m]
end