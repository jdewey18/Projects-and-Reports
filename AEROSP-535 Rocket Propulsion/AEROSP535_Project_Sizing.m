% AEROSP-535
% Author: Jacob Dewey
% Final Project Sizing Code
%--------------------------------------------------------------------------
% This code provides an estimate of the optimal stage and total mass of the
% vehicle using the design parameters (specific impulse and, structural 
% mass fractions, and intended delta-V). For any questions, please email
% deweyjm@umich.edu
%--------------------------------------------------------------------------

clc; clear;
g0 = 9.81; % [m/s^2] 
Isp_1 = 279.9; % Isp of First Stage Engines[s]
Isp_2 = 390;   % Isp of Second Stage Engines [s]
Isp_3 = 395;   % Isp of Third Stage Engine [s]
e1 = 0.05; % First Stage Mass Structure Fraction [kg/kg]
e2 = 0.07; % Second Stage Mass Structure Fraction [kg/kg]
e3 = 0.19; % Third Stage Mass Structure Fraction [kg/kg]
delta_v = 9500; % Targeted Delta-V for Mission [m/s]

% delta_v = (Isp_1 * g0 * log(R1)) + (Isp_2 * g0 * log(R2))

syms alpha
eqn1 = (Isp_1 * g0 * log(((alpha * Isp_1 * g0)+1)/(alpha * Isp_1 * g0 * e1))) + (Isp_2 * g0 * log(((alpha * Isp_2 * g0)+1)/(alpha * Isp_2 * g0 * e2))) + (Isp_3 * g0 * log(((alpha * Isp_3 * g0)+1)/(alpha * Isp_3 * g0 * e3))) == delta_v;
alpha = vpasolve(eqn1,alpha);

R1 = ((alpha * Isp_1 * g0)+1)/(alpha * Isp_1 * g0 * e1);
R2 = ((alpha * Isp_2 * g0)+1)/(alpha * Isp_2 * g0 * e2);
R3 = ((alpha * Isp_3 * g0)+1)/(alpha * Isp_3 * g0 * e3);

lambda1 = (1 - (e1*R1)) / (R1 - 1);
lambda2 = (1 - (e2*R2)) / (R2 - 1);
lambda3 = (1 - (e3*R3)) / (R3 - 1);

opt_mass_ratio = ((lambda1 + 1) / lambda1) * ((lambda2 + 1) / lambda2) * ((lambda3 + 1) / lambda3);

mL = 120000; % [kg] Payload Mass
m01 = opt_mass_ratio * mL; % [kg] Total Mass Including Payload

m_03 = double(mL * ((1 + lambda3) / lambda3));
m_02 = double(m_03 * ((1 + lambda2) / lambda2));
m_01 = double(m_02 * ((1 + lambda1) / lambda1));

m3 = m_03 - mL;
m2 = m_02 - m_03;
m1 = m_01 - m_02;

mb1 = double(m_01/R1);
mb2 = double(m_02/R2);
mb3 = double(m_03/R3);

g0 = 9.81;
T_W1 = 1.2;
T_W2 = 0.9;
T_W3 = 0.6;
% m_01 = 2616261.07;
% m_02 = 1554939.66;
% m_03 = 249557.23;

T1 = (T_W1 * g0 * m_01)/1000;
T2 = (T_W2 * g0 * m_02)/1000;
T3 = (T_W3 * g0 * m_03)/1000;

T_JMT1 = 6913; % [N] Thrust at sea-level
T_JMT2 = 1011.3; % [N] Thrust at vacuum
num_eng1 = T1 / T_JMT1;
num_eng2 = T2 / T_JMT2;
num_eng3 = T3 / T_JMT2;
T_act1 = 5 * T_JMT1;
T_act2 = 16 * T_JMT2;
T_act3 = 2 * T_JMT2;

Isp_JMT1 = 279.9; % [s] Isp at sea-level
Isp_JMT2 = 395.4; % [s] Isp at vacuum
m_dot_1act = T_act1 / (Isp_JMT1 * g0) * 1000;
m_dot_2act = T_act2 / (Isp_JMT2 * g0) * 1000;
m_dot_3act = T_act3 / (Isp_JMT2 * g0) * 1000;

% m_b1 = 1608005.7270020185076918395370261;
% m_b2 = 340933.99693052155398918675696194;
% m_b3 = 144615.87310254638997976755377669;
m_p1 = m_01 - mb1;
m_p2 = m_02 - mb2;
m_p3 = m_03 - mb3;
t_b1 = m_p1 / m_dot_1act;
t_b2 = m_p2 / m_dot_2act;
t_b3 = m_p3 / m_dot_3act;

fprintf('Stage 1 Wet Mass = %1.2f kg \n',m_01)
fprintf('Stage 1 Dry Mass = %1.2f kg \n',mb1)
fprintf('Stage 1 Thrust = %1.2f kN \n',T_act1)
fprintf('Stage 1 mdot = %1.2f kg/s \n',m_dot_1act)
fprintf('Stage 1 Burn Time = %1.2f s \n\n',t_b1)
fprintf('Stage 2 Wet Mass = %1.2f kg \n',m_02)
fprintf('Stage 2 Dry Mass = %1.2f kg \n',mb2)
fprintf('Stage 2 Thrust = %1.2f kN \n',T_act2)
fprintf('Stage 2 mdot = %1.2f kg/s \n',m_dot_2act)
fprintf('Stage 2 Burn Time = %1.2f s \n\n',t_b2)
fprintf('Stage 3 Wet Mass = %1.2f kg \n',m_03)
fprintf('Stage 3 Dry Mass = %1.2f kg \n',mb3)
fprintf('Stage 3 Thrust = %1.2f kN \n',T_act3)
fprintf('Stage 3 mdot = %1.2f kg/s \n',m_dot_3act)
fprintf('Stage 3 Burn Time = %1.2f s \n',t_b3)


