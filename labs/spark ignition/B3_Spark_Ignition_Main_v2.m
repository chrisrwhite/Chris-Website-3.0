%       Group: Thursday AM, B
%     Members: Kevin Myers, David Bedding, Samuel Hordeski,
%              Jorge Godoy, Justin Sandler, Chris White
%        Date: 3/22/15
%  Assignment: Spark Ignition Lab
% Description: Main Spark Ignition Program

%%  Clear memory and close windows
clear all
close all
clc

%% Sampling Rate
% Determined sampling rate of 250 or more points per second in order to capture the shape of the wave
% but actually collect at a sampling rate of the slowest piece (measuring
% device or board) 2300 samples per second

%%  Low Throttle
path = 'P:\me_475_3b\Spark Ignition\LowestThrottle\';
filename = 'SI_EngineData_VaryingLoadLT_Trial_';
startindex = 1;
endindex   = 55;
% sampling_time = 1/28000;

for j = startindex:endindex 
   data = dlmread([path filename num2str(j-1) '.xls']); 
   Crank_Angle_LT(:,j) = data(:,1);
   Voltage_LT(:,j) = data(:,2);
   Torque_LT(j) = mean(data(:,3));
   RPM_LT(j) = mean(data(:,4));
end

%   Volts --> psi
pressure_LT = 14.70 + 22.717*(Voltage_LT - 0.5);

%   Volume (based on cylinder motion)
a = 1.3;                          %crankshaft web length (in)
b = 4.134;                        %connecting rod length (in)
bore = 2.68;                      %(in)
clearance = 0.346;                %(in)
y_LT = Crank_Angle_LT;                               %crank angle
d_2_LT = a * cosd(y_LT);
c_LT = sqrt(a.^2 - d_2_LT.^2);
d_1_LT = sqrt(b.^2 - c_LT.^2);
d_3_LT = d_1_LT + d_2_LT;
volume_LT = (clearance + a + b - d_3_LT)*(pi*(bore^2)/4);

[b_,a_] = butter(3,.01);
pressure_LT_filt = filtfilt(b_,a_,pressure_LT);

%   Plot PV diagram
figure(1)
plot(volume_LT(:,2),pressure_LT_filt(:,2),'g')
title('Low Throttle')
xlabel('Volume (in^3)');
ylabel('Pressure (psi)');

%   Net work extracted from fluid
% work = zeros(1:5000);
% for n = 1:5000 
% work(n) = pressure(n) .* (volume(n+1) - volume(n));
% end
% power_fluid = diff(work)./(sampling_time);
% power_fluid_ft = power_fluid/12;
% power_hp_fluid = power_fluid_ft / 550.000037;

%   Low Throttle Shaft Power (Brake Work)
power_shaft_LT = Torque_LT .* RPM_LT;
power_shaft_LT_hp = power_shaft_LT / 5252;

%   Plot power vs speed
figure(4)
plot(RPM_LT,power_shaft_LT_hp,'g.')

%%  Medium Throttle
path = 'P:\me_475_3b\Spark Ignition\VaryingLoad (Mid Throttle)\';
filename = 'SI_EngineData_VaryingLoad_Trial_';
startindex = 1;
endindex   = 63;
% sampling_time = 1/28000;

for j = startindex:endindex 
   data = dlmread([path filename num2str(j-1) '.xls']); 
   Crank_Angle_MT(:,j) = data(:,1);
   Voltage_MT(:,j) = data(:,2);
   Torque_MT(j) = mean(data(:,3));
   RPM_MT(j) = mean(data(:,4));
end

%   Volts --> psi
Pressure_MT = 14.70 + 22.717*(Voltage_MT - 0.5);

%   Volume (based on cylinder motion)
y = Crank_Angle_MT;                  %crank angle
d_2 = a * cosd(y);
c = sqrt(a.^2 - d_2.^2);
d_1 = sqrt(b.^2 - c.^2);
d_3 = d_1 + d_2;
Volume_MT = (clearance + a + b - d_3)*(pi*(bore^2)/4);

Pressure_MT_filt = filtfilt(b_,a_,Pressure_MT);

%   Plot PV diagram
figure(2)
plot(Volume_MT(:,2),Pressure_MT_filt(:,2))
title('Med Throttle')
xlabel('Volume (in^3)');
ylabel('Pressure (psi)');

%   Net work extracted from fluid
% work = zeros(1:5000);
% for n = 1:5000 
% work(n) = pressure(n) .* (volume(n+1) - volume(n));
% end
% power_fluid = diff(work)./(sampling_time);
% power_fluid_ft = power_fluid/12;
% power_hp_fluid = power_fluid_ft / 550.000037;

%   Med Throttle Shaft Power (Brake Work)
Power_shaft_MT = Torque_MT .* RPM_MT;
Power_shaft_MT_hp = Power_shaft_MT / 5252;

%   Plot power vs speed
figure(4)
hold on
plot(RPM_MT,Power_shaft_MT_hp,'b.')

%%  Full Throttle
path = 'P:\me_475_3b\Spark Ignition\FullThrottle\';
filename = 'SI_EngineData_VaryingLoadFT_Trial_';
startindex = 1;
endindex   = 47;
% sampling_time = 1/28000;

for j = startindex:endindex 
   data = dlmread([path filename num2str(j-1) '.xls']); 
   Crank_Angle_FT(:,j) = data(300:end,1);
   Voltage_FT(:,j) = data(300:end,2);
   Torque_FT(j) = mean(data(:,3));
   RPM_FT(j) = mean(data(:,4));
end

%   Volts --> psi
pressure_FT = 14.70 + 22.717*(Voltage_FT - 0.5);

%   Volume (based on cylinder motion)
y_FT = Crank_Angle_FT;                               %crank angle
d_2_FT = a * cosd(y_FT);
c_FT = sqrt(a.^2 - d_2_FT.^2);
d_1_FT = sqrt(b.^2 - c_FT.^2);
d_3_FT = d_1_FT + d_2_FT;
volume_FT = (clearance + a + b - d_3_FT)*(pi*(bore^2)/4);

[b_,a_] = butter(3,.05);
pressure_FT_filt = filtfilt(b_,a_,pressure_FT);

%   Plot PV diagram
figure(3)
plot(volume_FT(:,2),pressure_FT_filt(:,2),'r')
title('Full Throttle')
xlabel('Volume (in^3)');
ylabel('Pressure (psi)');

%   Net work extracted from fluid
% work = zeros(1:5000);
% for n = 1:5000 
% work(n) = pressure(n) .* (volume(n+1) - volume(n));
% end
% power_fluid = diff(work)./(sampling_time);
% power_fluid_ft = power_fluid/12;
% power_hp_fluid = power_fluid_ft / 550.000037;

%   Full Throttle Shaft Power (Brake Work)
power_shaft_FT = Torque_FT .* RPM_FT;
power_shaft_FT_hp = power_shaft_FT / 5252;

%   Plot power vs speed
figure(4)
hold on
plot(RPM_FT,power_shaft_FT_hp,'r.')
title('Power vs Rotational Speed')
xlabel('RPM')
ylabel('Shaft Power (HP)')
legend('Low Throttle','Medium Throttle','Full Throttle','Location','Best')

%% Compare with given Briggs-Stratton Hp-RPM data (see excel sheet for formula)

BriggsStratton_RPM =[
1740
2000
2400
2800
3060
3200
3600
];

BriggsStratton_HP =[
7.05
8.43
10.35
11.95
12.74
13.07
13.71
];

% BriggsStratton_HP_Model = -1*10^(-6)*BriggsStratton_RPM.^2 + 0.0104*BriggsStratton_RPM - 7.3091;
% 
% BriggsStrattonComp(BriggsStratton_HP,BriggsStratton_HP_Model,BriggsStratton_RPM,RPM_FT) %Matrix fuckup here, i'm coming back




%% Deliverables:
% P-volume curves for selected operating conditions
% Relevant quantities vs. rpm for fixed throttle positions.
    % plot of efficiency of otto vs compression ratio?
% Mechanical efficiency calculated for selected operating conditions
% Comparison and discussion with Otto cycle idealization and provided data

%% Ideal Otto Cycle Comparison

% possible things to analyze:
% (a) the maximum temperature and pressure that occur during the cycle, 
% (b) the net work output, 
% (c) the thermal efficiency, and id) the mean effective pressure for the cycle.


%efficiency of otto cycle
Cp=1.005; %kJ/kg*K air @300 k
Cv=0.718; %kJ/kg*Kair @300 k
k=Cp/Cv;

%low throttle max and min volumes
V_BDC_LT= max(volume_LT(:,2));
V_TDC_LT=min(volume_LT(:,2));
r_LT= V_BDC_LT/V_TDC_LT
eff_otto_LT = 1 - 1/r_LT^(k-1)
%medium throttle max and min volumes
V_BDC_MT= max(Volume_MT(:,2));
V_TDC_MT=min(Volume_MT(:,2));
r_MT= V_BDC_MT/V_TDC_MT
eff_ott_MT = 1 - 1/r_MT^(k-1)
%high throttle max and min volumes
V_BDC_FT= max(volume_FT(:,2));
V_TDC_FT=min(volume_FT(:,2));
r_FT= V_BDC_FT/V_TDC_FT
eff_ott_FT = 1 - 1/r_FT^(k-1)

%  plot of efficiency of otto vs compression ratio THIS WAS POINTLESS ALL
%  THE SAME
figure(5)
plot(r_LT,eff_otto_LT,'r.',r_MT,eff_ott_MT,'bo',r_FT,eff_ott_FT,'gv')
title('efficiency of otto vs compression ratio')
xlabel('compression ratio')
ylabel('efficiency of otto')
legend('Low Throttle','Medium Throttle','Full Throttle','Location','Best')

%power and torque curves at MT vs speed
figure
plotyy(RPM_MT,Power_shaft_MT_hp,RPM_MT,Torque_MT)
title('Med Throttle')
xlabel('Speed (RPM)');
ylabel('Power (hp)');
legend('Medium Throttle Power','Medium Throttle Torque','Location','Best')

%conversions to SI
Pressure_MT=Pressure_MT*6.89475729;	 %psi -> kpa
Volume_MT(:,2)=Volume_MT(:,2)*1.63871e-5; %in^3 -> m^3


%% thermodynamic analysis of air standard otto cycle page 87 of pdf
% 6-1-2-3-4-5-6 (6 being TDC intake)

%ideal gas law find temperatures throughout cycle
R = 8.314; %J/mol K - ideal gas constant 
temp_MT=Pressure_MT(:,2).*Volume_MT(:,2)/R; 

%process 6-1: isobaric intake of air @ P_o, intake valve open and exaust
%valve closed
P_o = min(Pressure_MT(:,2));
P_6 = P_o;
P_1 = P_6; %maybe???????????

V_BDC_MT= max(Volume_MT(:,2));
V_TDC_MT=min(Volume_MT(:,2));
v_1=V_BDC_MT; %m^3/kg - specific volume
v_6 =V_TDC_MT; %m^3/kg - specific volume

w_61 = P_o*(v_1 - v_6);
v_2 = v_6;
%process 1-2: isentropic compression stroke @ P_1, valves closed

T_1 =P_1*v_1/R;
r_c = v_1/v_2; %compresson ratio

T_2 = T_1*(r_c)^(k-1);
P_2 = P_1*(r_c)^k;

w_12 = Cv*(T_1 - T_2);

%process 2-3: isobaric combustion stroke, valves closed
v_TDC = min(Volume_MT(:,2));

v_3= v_2;
v_6 =v_3;

T_3 = max(temp_MT); %define temp array
P_3 = max(Pressure_MT(:,2));
q_23 = Cv*(T_3-T_2); %q in


%process 3-4: isentropic  power stroke, valves closed
T_4 = T_3*(1/r_c)^(k-1);
P_4 = P_3*(1/r_c)^k;
w_34 = Cv*(T_3 - T_4);

%process 4-5: constant volume heat rejection  exhaust valve open intake
%closed
v_BDC = max(Volume_MT(:,2));
v_1 =v_BDC;
v_4 =v_1;
v_5 =v_4;
q_45 =Cv*(T_1-T_4);%q out

%process 5-6: isobaric exaust stroke with  exhaust valve open intake
%closed
P_5= P_6;
w_56 = P_o*(v_6 - v_5);

%efficiency of otto cycle
eff_otto = 1 - (1/r_c)^(k-1);






