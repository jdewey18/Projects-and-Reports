% AEROSP-530
% Jacob Dewey
% Final Project
%--------------------------------------------------------------------------
% This code performs a cycle analysis for the decided configuration of the
% turbojet engine that meets the requirements of the RFP. Much of the data
% below in regards to the efficiencies, etc is empirical data of previous
% engines. The file is self-contained and has all of the required functions
% for its execution. 
% 
% Please direct any questions to deweyjm@umich.edu
%--------------------------------------------------------------------------

clc, clear all, close all


%% Problem Constants:
% Flight Conditions and Fluid Properties
amb.Alt = 30e3; % Altitude [ft]
amb.tBurn = 110*60; % Time at Cruise [s];
amb.fuelDens = 6.532; % Fuel Density [lbm/gal]
amb.Ma = 0.8; % Freestream Mach Number
amb.M1 = 0.4; % Mach Number at Inlet Exit
amb.Pa = 4.36; % Ambient Pressure [psia]
amb.Ta = convtemp(-44.4,'C','R'); % Ambient Temperature [Rankine]
amb.gamc = 1.4; % Specific Heat Ratio (cold)
amb.gamh = 1.3333; % Specific Heat Ratio (hot)
amb.R = 53.35; % Gas Constant [(ft-lbf)/(lbm-R)]
amb.J = 778; % Energy Constant [(ft-lbf)/Btu]
amb.gc = 32.174; % Gravitational Constant [(ft-lbm)/(lbf-s^2)]
amb.LHV = 18550; % Fuel Lower Heating Value [Btu/(lbm-fuel)]
amb.Fsg = 0.8165; % Fuel Specific Gravity
amb.Cpc = amb.R/(1-1/amb.gamc)/amb.J; % Specific Heat (cold) [(ft-lbf)/(lbm-R)]
amb.Cph = amb.R/(1-1/amb.gamh)/amb.J; % Specific Heat (hot) [(ft-lbf)/(lbm-R)]

% Component Information:
comp.nr = 0.99; % Inlet Recovery (Pt1/Pt0)
comp.ncf = 0.88; % Fan Efficiency
comp.nc = 0.87; % Compressor Efficiency
comp.mlm2 = 0.01; % Compressor Leakage Flow (m_leak/m_2)
comp.mcm2 = 0.06; % Turbine Cooling Flow (m_cool/m_2)
comp.mbm2 = 0.03; % Customer Bleed Flow (m_bleed/m_2)
comp.Pt4Pt3 = 0.95; % Combustor Pressure Drop (Pt4/Pt3.1)
comp.nb = 0.995; % Combustor Efficiency
comp.nt = 0.92; % Core Turbine Efficiency
comp.ntf = 0.925; % Fan Turbine Efficiency
comp.Cv = 0.983; % Core Exhaust Velocity Coefficient
comp.Cvf = 0.985; % Fan Exhaust Velocity Coefficient

% Limitations and Other:
lim.maxTt3 = convtemp(450,'F','R'); % Max Compressor Discharge Temperature [R]
lim.maxTt495 = convtemp(890,'C','R'); % Max Exhaust Gas Temperature (EGT) [R]
lim.maxTt9 = convtemp(880,'C','R'); % Max Core Exhaust Nozzle Gas Temperature [R]
lim.maxD1 = 30; % Maximum Internal Flowpath Diameter [in]
lim.maxFuelVol = 1600; % Maximum Fuel [gal]
lim.maxFuelm = amb.fuelDens*lim.maxFuelVol;

%% Cycle Analysis:
var.BPR = 0; % Bypass Ratio (m3/m2)
var.FPR = 1; % Fan Pressure Ratio (Pt13/Pt1 and Pt2/Pt1)
var.CPR = 8; % Compressor Pressure Ratio (Pt3/Pt2)
var.Tt4 = convtemp(2100,'F','R'); % Combustor Exit Temperature [Rankine]
var.target = 4800;
var.text = "Design: Configuration 1";
[Config1,var] = TurbojetConfig(amb,comp,var);
disp(' ')
dispSolution(Config1,var.BPR,var.text);


%% Functions:
% Turbofan/Turbojet Cycle Analysis
function [out] = TurbofanSim(amb,comp,var)
% Station 0 (Inlet Plane):
out.Pt0 = amb.Pa*(1 + (amb.gamc-1)*amb.Ma^2/2)^(amb.gamc/(amb.gamc-1)); % Total Pressure at Station 0
out.Tt0 = amb.Ta*(1 + (amb.gamc-1)*amb.Ma^2/2); % Total Temperature at Station 0
out.m0 = var.mdota; % Mass Flow Rate at Station 0

% Station 1 (Diffuser Exit):
out.Pt1 = out.Pt0*comp.nr;
out.Tt1 = out.Tt0;
out.m1 = out.m0;
out.Ts1 = out.Tt1/(1 + (amb.gamc - 1)*amb.M1^2/2);
out.Ps1 = out.Pt1/(1 + (amb.gamc - 1)*amb.M1^2/2)^(amb.gamc/(amb.gamc-1));
out.a1 = sqrt(amb.gamc*amb.gc*amb.R*out.Ts1);
out.u1 = amb.M1*out.a1;
out.rhos1 = out.Ps1*144/(amb.R*out.Ts1);
out.A1 = out.m1/(out.rhos1*out.u1);

% Station 2 (Fan Exit):
out.Pt2 = out.Pt1*var.FPR;
out.Tt2 = out.Tt1*((1/comp.ncf)*(var.FPR^((amb.gamc-1)/amb.gamc) - 1) + 1);
out.m2 = out.m1/(var.BPR + 1);

% Station 3 (Compressor Exit Before Bleeds):
out.Pt3 = out.Pt2*var.CPR;
out.Tt3 = out.Tt2*((1/comp.nc)*(var.CPR^((amb.gamc-1)/amb.gamc) - 1) + 1);
out.m3 = out.m2;

% Station 3.1 (Compressor Exit After Bleeds):
out.Pt31 = out.Pt3;
out.Tt31 = out.Tt3;
out.mcool = out.m2*comp.mcm2;
out.m31 = out.m3*(1 - comp.mbm2 - comp.mlm2) - out.mcool;

% Station 4 (Combustor Exit)
out.Pt4 = out.Pt3*comp.Pt4Pt3;
out.Tt4 = var.Tt4;
out.mf = (out.m31*amb.Cph*(out.Tt4-out.Tt31))/(amb.LHV*comp.nb); %-amb.Cph*out.Tt4);
out.m4 = out.m31 + out.mf;

% Station 4.9 (High-Pressure Turbine Exit):
out.m49 = out.m4;
out.Tt49 = out.Tt4 - out.m2*amb.Cpc*(out.Tt3 - out.Tt2)/(out.m4*amb.Cph);
% out.Tt49 = out.Tt4 + (out.Tt49i - out.Tt4)*comp.nt;
out.Pt49 = out.Pt4*(1 - (1 - out.Tt49/out.Tt4)/comp.nt)^(amb.gamh/(amb.gamh - 1));

% Station 4.95 (Low Pressure Turbine Entrance):
out.m495 = out.m49 + out.mcool;
out.Tt495 = (out.m49*amb.Cph*out.Tt49 + out.mcool*amb.Cpc*out.Tt3)/(out.m495*amb.Cph);
out.Pt495 = out.Pt49; 

% Station 5 (Low Pressure Turbine Exit): 
out.m5 = out.m495;
out.Tt5 = out.Tt495 - (out.m1*amb.Cpc*(out.Tt2-out.Tt1))/(out.m5*amb.Cph);
% out.Tt5 = out.Tt495 + (out.Tt5i - out.Tt495)*comp.ntf;
out.Pt5 = out.Pt495*(1 - (1 - out.Tt5/out.Tt495)/comp.ntf)^(amb.gamh/(amb.gamh - 1));

% Station 9; (Nozzle Exit):
out.m9 = out.m5;
out.Tt9 = out.Tt5;
out.P9 = amb.Pa;
out.T9i = out.Tt9*(out.P9/out.Pt5)^((amb.gamh-1)/amb.gamh);
out.T9 = out.Tt5 - comp.Cv^2*(out.Tt5-out.T9i);
out.u9 = sqrt(2*amb.Cph*amb.gc*amb.J*(out.Tt5 - out.T9));
out.M9 = out.u9/sqrt(amb.gamh*amb.gc*amb.R*out.T9);
out.Pt9 = out.P9*(1 + 0.5*(amb.gamh-1)*out.M9^2)^(amb.gamh/(amb.gamh-1));

% Station 13 (Bypass Fan Exit):
out.Pt13 = out.Pt1*var.FPR;
out.Tt13i = out.Tt1*(var.FPR^((amb.gamc-1)/amb.gamc));
out.Tt13 = out.Tt1 + (out.Tt13i - out.Tt1)/comp.ncf;
out.m13 = var.BPR*out.m1/(var.BPR + 1);

% Station 19 (Fan Nozzle Exit):
out.m19 = out.m13;
out.Tt19 = out.Tt13;
out.P19 = amb.Pa;
out.T19i = out.Tt19*(out.P19/out.Pt13)^((amb.gamc-1)/amb.gamc);
out.T19 = out.Tt13 - comp.Cvf^2*(out.Tt13-out.T19i);
out.u19 = sqrt(2*amb.Cpc*amb.gc*amb.J*(out.Tt13 - out.T19));
out.M19 = out.u19/sqrt(amb.gamc*amb.gc*amb.R*out.T19);
out.Pt19 = out.P19*(1 + 0.5*(amb.gamc-1)*out.M19^2)^(amb.gamc/(amb.gamc-1));

% Processing:
out.D1 = 2*12*sqrt(out.A1/pi);
out.u0 = amb.Ma*sqrt(amb.R*amb.gamc*amb.gc*amb.Ta);
out.thrust = (out.m9*out.u9 + out.m19*out.u19 - out.m0*out.u0)/amb.gc;
out.specThrust = out.thrust/out.m2;
out.fuelBurned = out.mf*amb.tBurn;
end

% Finds the Mass Flow That Satisfies the Target Thrust Value
function [out,var] = TurbojetConfig(amb,comp,var)
mdotvec = [0,200];
margin = 1e-2;
out.thrust = 0;
while abs(out.thrust-var.target) > margin
    var.mdot2 = mean(mdotvec);
    var.mdota = var.mdot2*(var.BPR+1);
    [out] = TurbofanSim(amb,comp,var);
    
    if out.thrust > var.target+margin
        mdotvec(2) = var.mdot2;
    elseif out.thrust < var.target - margin
        mdotvec(1) = var.mdot2;
    end
end
end

% Print Outputs:
function dispSolution(out,BPR,text)
disp('-----------------------------------')
disp(text)
disp('-----------------------------------')
disp('           Tt      pt       mdot')
disp('Station    (R)    (psia)   (lbm/s)')
disp('__________________________________')
fprintf('0        %4.2f    %4.2f     %4.2f \n',out.Tt0,out.Pt0,out.m0)
fprintf('1        %4.2f    %4.2f     %4.2f \n',out.Tt1,out.Pt1,out.m1)
fprintf('2        %4.2f    %4.2f     %4.2f \n',out.Tt2,out.Pt2,out.m2)
fprintf('3        %4.2f    %4.2f     %4.2f \n',out.Tt3,out.Pt3,out.m3)
fprintf('3.1      %4.2f    %4.2f     %4.2f \n',out.Tt31,out.Pt31,out.m31)
fprintf('4        %4.2f    %4.2f     %4.2f \n',out.Tt4,out.Pt4,out.m4)
fprintf('4.9      %4.2f    %4.2f     %4.2f \n',out.Tt49,out.Pt49,out.m49)
fprintf('4.95     %4.2f    %4.2f     %4.2f \n',out.Tt495,out.Pt495,out.m495)
if BPR > 0
    fprintf('5        %4.2f    %4.2f     %4.2f \n',out.Tt5,out.Pt5,out.m5)
end
fprintf('9        %4.2f    %4.2f     %4.2f \n',out.Tt9,out.Pt9,out.m9)

if BPR > 0
    fprintf('13       %4.2f    %4.2f     %4.2f \n',out.Tt13,out.Pt13,out.m13)
    fprintf('19       %4.2f    %4.2f     %4.2f \n',out.Tt19,out.Pt19,out.m19)
end
disp(' ')
fprintf('u9 =  %1.2f ft/s\n',out.u9)
if BPR > 0
fprintf('u19 = %1.2f ft/s\n',out.u19)
end
fprintf('Thrust: %1.2f lbf \n',out.thrust)
fprintf('Specific Thrust: %1.2f lbf-s/lbm \n',out.specThrust)
fprintf('Fuel Burned: %1.2f lbm \n',out.fuelBurned)
fprintf('Max Diameter: %1.2f in \n',out.D1)
fprintf('Tt3:   %1.2f R \n',out.Tt3)
fprintf('Tt495: %1.2f R \n',out.Tt495)
fprintf('Tt9:   %1.2f R \n',out.Tt9)
end

