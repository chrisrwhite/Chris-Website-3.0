%       Group: Thursday AM, B
%     Members: Kevin Myers, David Bedding, Samuel Hordeski,
%              Jorge Godoy, Justin Sandler, Chris White
%        Date: 2/8/15
%  Assignment: Wind Tunnel 1

%%  Clear memory and close windows
clear all
close all
clc
%%  Read published Schlichting sphere data
SphereData_Sch = xlsread('Schlichting Sphere Data.xlsx');
Re_Sch = SphereData_Sch(:,1);
Cd_Sch = SphereData_Sch(:,2);

%%  Define Known Values
%   Sphere radii
D_15 = 1.673*.0254;
D_3 = 2.935*.0254;
D_4 = 4*.0254;
D_Golf = 1.69*.0254;

%   Force Readout Scaling Factor
Calibration_Factor = 95.5;

%   Fluid density
rho = 1.20;

%   Dynamic Viscosity
mu = 1.825*10^(-5);

%%  Read experimantal data
%   Sphere, D = 1.5in
SphereData_15 = xlsread('Sphere Data D=1.5.xlsx');
P_air_15 = SphereData_15(:,1);
P_air_dyn_15 = -P_air_15;
F_drag_uncalibrated_15 = SphereData_15(:,2);
F_drag_15 = F_drag_uncalibrated_15 / Calibration_Factor;
V_air_15 = sqrt(P_air_dyn_15 * 2 / rho);

%   Sphere, D = 3in
SphereData_3 = xlsread('Sphere Data D=3.xlsx');
P_air_3 = SphereData_3(:,1);
P_air_dyn_3 = -P_air_3;
F_drag_uncalibrated_3 = SphereData_3(:,2);
F_drag_3 = F_drag_uncalibrated_3 / Calibration_Factor;
V_air_3 = sqrt(P_air_dyn_3 * 2 / rho);

%  Sphere, D = 4in
SphereData_4 = xlsread('Sphere Data D=4.xlsx');
P_air_4 = SphereData_4(:,1);
P_air_dyn_4 = -P_air_4;
F_drag_uncalibrated_4 = SphereData_4(:,2);
F_drag_4 = F_drag_uncalibrated_4 / Calibration_Factor;
V_air_4 = sqrt(P_air_dyn_4 * 2 / rho);

%   Sphere, D = 4in Turbulent Trip
SphereData_4_turb = xlsread('Sphere Data D=4 TurbTrip.xlsx');
P_air_4_turb = SphereData_4_turb(:,1);
P_air_dyn_4_turb = -P_air_4_turb;
F_drag_uncalibrated_4_turb = SphereData_4_turb(:,2);
F_drag_4_turb = F_drag_uncalibrated_4_turb / Calibration_Factor;
V_air_4_turb = sqrt(P_air_dyn_4_turb * 2 / rho);

%   Golf Ball
SphereData_Golf = xlsread('Sphere Data Golfball.xlsx');
P_air_Golf = SphereData_Golf(:,1);
P_air_dyn_Golf = -P_air_Golf;
F_drag_uncalibrated_Golf = SphereData_Golf(:,2);
F_drag_Golf = F_drag_uncalibrated_Golf / Calibration_Factor;
V_air_Golf = sqrt(P_air_dyn_Golf * 2 / rho);

%%  Calculate Reynolds Numbers
Re_15 = (rho*D_15*V_air_15)/mu;
Re_3 = (rho*D_3*V_air_3)/mu;
Re_4 = (rho*D_4*V_air_4)/mu;
Re_4_turb = (rho*D_4*V_air_4_turb)/mu;
Re_Golf = (rho*D_Golf*V_air_Golf)/mu;

%%  Calculate Coefficiants of Drag
%   Projected frontal areas
A_15 = (pi/4)*D_15^2;
A_3 = (pi/4)*D_3^2;
A_4 = (pi/4)*D_4^2;
A_Golf = (pi/4)*D_Golf^2;

%   Coefficients of drag
Cd_15 = F_drag_15 ./(P_air_dyn_15*A_15);
Cd_3 = F_drag_3 ./(P_air_dyn_3*A_3);
Cd_4 = F_drag_4 ./(P_air_dyn_4*A_4);
Cd_4_turb = F_drag_4_turb ./(P_air_dyn_4_turb*A_4);
Cd_Golf = F_drag_Golf ./(P_air_dyn_Golf*A_Golf);

%   Coefficients of drag modeled for Schlichting sphere
Cd_model_Sch = 24 ./ Re_Sch + (2.6 .* (Re_Sch ./ 5.0)) ./(1 + (Re_Sch ./ 5.0) .^ 1.52) + (0.411 .* (Re_Sch ./ 263000) .^ (-7.94)) ./ (1+(Re_Sch ./ 263000) .^(-8)) + (Re_Sch .^(.80)) ./ 461000;

%%  Plot data
figure(1)
loglog(Re_Sch,Cd_Sch,'ok')
hold on
loglog(Re_15,Cd_15,'xg')
loglog(Re_3,Cd_3,'*r')
loglog(Re_4,Cd_4,'^b')
legend('Schlichting','1.5', '3','4','Location','Southeast')
xlim([10^4 5*10^5])
title('Smooth Spheres')

figure(2)
loglog(Re_Sch,Cd_Sch,'-ok')
hold on
loglog(Re_15,Cd_15,'-xg')
loglog(Re_3,Cd_3,'-*r')
loglog(Re_4,Cd_4,'-^b')
loglog(Re_4_turb,Cd_4_turb,'-sm')
loglog(Re_Golf,Cd_Golf,'-+k')
legend('Schlichting', '1.5" Smooth', '3" Smooth','4" Smooth','4" Tripped','Golf Ball','Location','Southeast')
xlim([10^4 5*10^5])
ylim([12^-1 10^0])
xlabel('Reynolds Number')
ylabel('Coefficient of Drag')
title('Smooth Spheres & Tripped Spheres')

figure(3)
loglog(Re_15,Cd_15,'.b');
hold on
loglog(Re_3,Cd_3,'.r');
loglog(Re_4,Cd_4,'.g');
loglog(Re_Sch,Cd_Sch,'mx')
loglog(Re_Sch,Cd_model_Sch,'k');
xlim([10^4 5*10^5])
ylim([10^-.75 10^-.1])
xlabel('Reynolds Number')
ylabel('Coefficient of Drag')
legend('1.5" Sphere','3" Sphere','4" Sphere','Schlichting','Morrison Correlation Model','Location','Best')
title('Calculated Cd vs. Modeled Cd for Smooth Spheres')

figure(4)
loglog(Re_Sch,Cd_Sch, '.k',Re_Sch,Cd_model_Sch,'r')
xlabel('Reynolds Number')
ylabel('Coefficient of Drag')
title('Given Schlichting Data vs. Morrison Correlation')

%%  Read Provided Airfoil Data "NACA Results"
Airfoil_Data_Provided = xlsread('Provided Airfoil Data.xlsx','NACA Results');
Alpha_Prov = Airfoil_Data_Provided(:,1);
Cl_Prov = Airfoil_Data_Provided(:,2);
Cd_Prov = Airfoil_Data_Provided(:,3);

%%  Read Provided Airfoil Data "Computational Results"
Xfoil_Data_Provided = xlsread('Provided Airfoil Data.xlsx',1);
Alpha_Xfoil_Re1e5 = Xfoil_Data_Provided(:,1);
Cl_Xfoil_Re1e5 = Xfoil_Data_Provided(:,2);
Cd_Xfoil_Re1e5 = Xfoil_Data_Provided(:,3);
Alpha_Xfoil_Re5e5 = Xfoil_Data_Provided(:,5);
Cl_Xfoil_Re5e5 = Xfoil_Data_Provided(:,6);
Cd_Xfoil_Re5e5 = Xfoil_Data_Provided(:,7);

%%  Plot Provided Airfoil Data
figure(5)
plot(Alpha_Prov,Cl_Prov,'+-k', Alpha_Prov,Cd_Prov,'x-k')
xlabel('Angle of Attack (deg)')
legend('Coefficient of Lift (NACA)','Coefficient of Drag (NACA)','Location','Best')
title('NACA Airfoil Data (Re = 8e5)')

%%  
figure(6)
hold on
plot(Alpha_Xfoil_Re1e5,Cl_Xfoil_Re1e5,'b', Alpha_Xfoil_Re1e5,Cd_Xfoil_Re1e5,'m')
plot(Alpha_Xfoil_Re5e5,Cl_Xfoil_Re5e5,'-.b', Alpha_Xfoil_Re5e5,Cd_Xfoil_Re5e5,'-.m')
xlabel('Angle of Attack (deg)')
legend('Coefficient of Lift, Re = 1e5','Coefficient of Drag Re = 1e5','Coefficient of Lift, Re = 5e5','Coefficient of Drag Re = 5e5','Location','Best')
title('Xfoil Data (Re = 1e5 & Re = 5e5)')

%% Read Experimental Airfoil Data
Airfoil_Data_Collected = xlsread('Airfoil Data.xlsx');
Norm_Collected_No_Vel = Airfoil_Data_Collected(:,1);        %Units:Counts
Axial_Collected_No_Vel = Airfoil_Data_Collected(:,2);       %Units:Counts
Alpha = Airfoil_Data_Collected(:,3);                        %Units:Degrees
Pressure_Collected_Top_Vel = Airfoil_Data_Collected(:,4);   %Units:Pascals
Norm_Collected_Top_Vel = Airfoil_Data_Collected(:,5);       %Units:Counts
Axial_Collected_Top_Vel = Airfoil_Data_Collected(:,6);      %Units:Counts

%   Other Info
V_airfoil = 39.5; %(m/s)
b=10*.0254; %(m)
c=3.498*.0254; %(m)
Planform_A = b*c; %(m^2)
Re_airfoil = (V_airfoil*rho*c)/mu;

%%  Corrected Normal & Axial Forces (for gravity)
%   Converting to Newtons
F_N_Top_Vel = Norm_Collected_Top_Vel/Calibration_Factor;
F_N_No_Vel = Norm_Collected_No_Vel/Calibration_Factor;
F_A_Top_Vel = Axial_Collected_Top_Vel/Calibration_Factor;
F_A_No_Vel = Axial_Collected_No_Vel/Calibration_Factor;

%   Corrected Forces
F_N_Correct = F_N_Top_Vel - F_N_No_Vel;
F_A_Correct = F_A_Top_Vel - F_A_No_Vel;

%%  Calculate Lift and Drag Forces
%   Convert from degrees to radians
Alpha_rad = Alpha*(pi/180);

F_L = F_N_Correct.*cos(Alpha_rad) - F_A_Correct.*sin(Alpha_rad);
F_D = F_N_Correct.*sin(Alpha_rad) + F_A_Correct.*cos(Alpha_rad);

%%  Calculate Coefficients of Lift & Drag
Cl_airfoil = F_L ./(Pressure_Collected_Top_Vel*Planform_A);
Cd_airfoil = F_D ./(Pressure_Collected_Top_Vel*Planform_A);

%%  Plot Airfoil Data
figure(7)
plot(Alpha,Cl_airfoil,'-+g', Alpha,Cd_airfoil,'-xr')
xlabel('Angle of Attack (deg)')
legend('Coefficient of Lift','Coefficient of Drag','Location','Best')
title('Airfoil Data, 2D Assumpion')

%%  Corrected Cl & Cd models for 3D airfoil
e_1 = 0.7; %rectangular wing planform

%   Aspect Ratio
AR = b/c;

Cl_airfoil_2D = 2*pi*(Alpha_rad-Alpha_rad(1));
Cl_airfoil_3D = Cl_airfoil_2D*AR/(AR+2/e_1);

Cd_0 = 0.025; %aerospaceweb.org
Cd_airfoil_3D = Cd_0+Cl_airfoil_2D.^2/(4*pi*AR*e_1);

%%  Plot Airfoil Data, Provided Data, and Prandtl Lifting Line Models
figure(8)
plot(Alpha,Cl_airfoil,'-+g', Alpha,Cl_airfoil_3D,'-b', Alpha_Prov,Cl_Prov,'-ok')
xlabel('Angle of Attack (deg)')
legend('Coefficient of Lift (Experimental)','Coefficient of Lift (3D Model)','Coefficient of Lift (NACA Data)','Location','Best')
title('Airfoil Lift Performance Comparison')

figure(9)
plot(Alpha,Cd_airfoil,'-xr', Alpha,Cd_airfoil_3D,'-m', Alpha_Prov,Cd_Prov,'-ok')
xlabel('Angle of Attack (deg)')
legend('Coefficient of Drag (Experimental)','Coefficient of Drag (3D Model)','Coefficient of Drag (NACA Data)','Location','Best')
title('Airfoil Drag Performance Comparison')

figure(10)
plot(Alpha,Cd_airfoil,'+', Alpha,Cd_airfoil_3D,'-o',Alpha,Cl_airfoil,'-*', Alpha,Cl_airfoil_3D,'-.',Alpha_Prov,Cl_Prov,'-x',Alpha_Prov,Cd_Prov,'s',Alpha_Xfoil_Re1e5,Cl_Xfoil_Re1e5,'d', Alpha_Xfoil_Re1e5,Cd_Xfoil_Re1e5,'^',Alpha_Xfoil_Re5e5,Cl_Xfoil_Re5e5,'v',Alpha_Xfoil_Re5e5,Cd_Xfoil_Re5e5,'h')
xlabel('Angle of Attack (deg)')
ylabel('Coefficient of Drag')
legend('COF of Drag (Experimental)','COF of Drag (3D Model)','COF of Lift (Experimental)','COF of Lift (3D Model)','COF of Lift (NACA Data)','COF of Drag (NACA Data)','COF of Lift, Re = 1e5 (Xfoil)','COF of Drag, Re = 1e5 (Xfoil)','COF of Lift, Re = 5e5 (Xfoil)','COF of Drag, Re = 5e5 (Xfoil)','Location','Best')
title('Airfoil Performance, Experimental vs Theoretical Lift & Drag')

figure(11)
plot(Alpha,Cl_airfoil,'-+g', Alpha,Cd_airfoil,'-xr',Alpha_Xfoil_Re1e5,Cl_Xfoil_Re1e5,'b', Alpha_Xfoil_Re1e5,Cd_Xfoil_Re1e5,'m',Alpha_Xfoil_Re5e5,Cl_Xfoil_Re5e5,'-.b', Alpha_Xfoil_Re5e5,Cd_Xfoil_Re5e5,'-.m')
xlabel('Angle of Attack (deg)')
legend('Coefficient of Lift (Experimental)','Coefficient of Drag (Experimental)','Coefficient of Lift, Re = 1e5','Coefficient of Drag, Re = 1e5','Coefficient of Lift, Re = 5e5','Coefficient of Drag, Re = 5e5','Location','Best')
title('Airfoil Performance, Experimental vs Xfoil Data (Re = 1e5 & Re = 5e5')
