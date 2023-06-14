function [T, P, rho] = atmo_model_100km(maxAlt)


% Defining constants
g0 = 9.81; % m/s^2, gravity at sea level
T0 = 288.16; % K, temp at sea level
P0 = 101325; % Pa, pressure at sea level
rE = 6378.154; % km, radius of Earth
R = 287.05; % J/kg-K, gas constant for air
rho0 = P0/(R*T0); % kg/m^3, density at sea level

% Defining lapse rates for gradient regions
a1 = -6.5e-3; % K/m
a2 = 3e-3; % K/m
a3 = -4.5e-3; % K/m
a4 = 4e-3; % K/m

for alt = 0:1:maxAlt
    % Gradient Region: 0 km to 11 km
    if 0 <= alt && alt <= 11
        T(alt+1) = T0 + (a1*alt*1000); % K
        P(alt+1) = P0*((T(alt+1)/T0)^((-g0/(a1*R)))); % Pa
        rho(alt+1) = rho0*((T(alt+1)/T0)^(-((g0/(a1*R))+1))); % kg/m^3
        
        % Isothermal Region: 11 km to 25 km
    elseif 12 <= alt && alt <= 25
        T(alt+1) = T(12); % K
        P(alt+1) = P(12)*exp((-g0/(R*T(alt+1)))*(alt-11)*1000); % Pa
        rho(alt+1) = rho(12)*exp((-g0/(R*T(alt+1)))*(alt-11)*1000); % kg/m^3
        
        % Gradient Region: 25 km to 47 km
    elseif 26 <= alt && alt <= 47
        T(alt+1) = T(26) + (a2*(alt-25)*1000); % K
        P(alt+1) = P(26)*((T(alt+1)/T(26))^-(g0/(a2*R))); % Pa
        rho(alt+1) = rho(26)*((T(alt+1)/T(26))^(-((g0/(a2*R))+1))); % kg/m^3
        
        % Isothermal Region: 47 km to 53 km
    elseif 48 <= alt && alt <= 53
        T(alt+1) = T(48); % K
        P(alt+1) = P(48)*exp((-g0/(R*T(alt+1)))*(alt-47)*1000); % Pa
        rho(alt+1) = rho(48)*exp((-g0/(R*T(alt+1)))*(alt-47)*1000); % kg/m^3
        
        % Gradient Region: 53 km to 79 km
    elseif 54 <= alt && alt <= 79
        T(alt+1) = T(54) + (a3*(alt-53)*1000); % K
        P(alt+1) = P(54)*((T(alt+1)/T(54))^-(g0/(a3*R))); % Pa
        rho(alt+1) = rho(54)*((T(alt+1)/T(54))^(-((g0/(a3*R))+1))); % kg/m^3
        
        % Isothermal Region: 79 km to 90 km
    elseif 80 <= alt && alt <= 90
        T(alt+1) = T(80); % K
        P(alt+1) = P(80)*exp((-g0/(R*T(alt+1)))*(alt-79)*1000); % Pa
        rho(alt+1) = rho(80)*exp((-g0/(R*T(alt+1)))*(alt-79)*1000); % kg/m^3
        
        % Gradient Region: 90 km to 100 km
    elseif 91 <= alt && alt <= 100
        T(alt+1) = T(91) + (a4*(alt-90)*1000); % K
        P(alt+1) = P(91)*((T(alt+1)/T(91))^-(g0/(a4*R))); % Pa
        rho(alt+1) = rho(91)*((T(alt+1)/T(91))^(-((g0/(a4*R))+1))); % kg/m^3
        
    end
end
% Dynamic Viscosity
mu(alt+1) = ((1.458e-6)*sqrt(T(alt+1)))/(1+(110.4/T(alt+1))); % Pa-s
end
