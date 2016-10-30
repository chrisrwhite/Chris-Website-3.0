%%   Filename: PreLab10.m
%  Assignment: Problem Set 5
%      Author: Michael Duncan
%        Date: 10/30/14
% Description: Graphs of root loci and step responses of 
% transfer functions from Lab 9

%close any open graphical windows and clear the variables
close all
clc
clear

s=tf('s');

%% Open Loop Transfer Functions
 G_1 = 8* (s+7.6)/(s^2+7.6*s+444);
 G_2 = 8* (s+11.4)*(s+7.6)/(s^2+7.6*s+444);
 G_3 = 7*1/14*-40/(s^2+0.44*s+80);
 G_4 = 7*1/14*-40*(s+18)/(s^2+0.44*s+80);
 
% G_msl1 = TransFun(1,3,1);
%% Root Loci
figure;
step(G_1)
hold on
step(G_2)

hold off
figure;
step(G_3)
hold on
step(G_4)
hold off

% figure(1)
% rlocus(G_m1)
% title('Root Locus (Motor Alone)')
% figure(2)
% rlocus(G_ml1)
% title('Root Locus (Motor and Load Inertia)')
% figure(3)
% rlocus(G_msl1)
% title('Root Locus (Motor, Spring and Load Inertia)')
% %% Responses
% % Kp = 1
% [R_m1,t_m1] = step(TransFun(1,1,0));
% [R_ml1,t_ml1] = step(TransFun(1,2,0));
% [R_msl1,t_msl1] = step(TransFun(1,3,0));
% % Kp = 2
% [R_m2,t_m2] = step(TransFun(2,1,0));
% [R_ml2,t_ml2] = step(TransFun(2,2,0));
% [R_msl2,t_msl2] = step(TransFun(2,3,0));
% % Kp = 5
% [R_m5,t_m5] = step(TransFun(5,1,0));
% [R_ml5,t_ml5] = step(TransFun(5,2,0));
% [R_msl5,t_msl5] = step(TransFun(5,3,0));
% %% Plot Responses
% figure(4)
% plot(t_m1,R_m1,'r-',t_m2,R_m2,'b--',t_m5,R_m5,'g-.')
% title('Step Responses with Varying Gain (Motor Alone)')
% xlabel('Time(s)')
% ylabel('Angular Velocity (rad/s)')
% legend('Kp=1','Kp=2','Kp=5')
% figure(5)
% plot(t_ml1,R_ml1,'r-',t_ml2,R_ml2,'b--',t_ml5,R_ml5,'g-.')
% title('Step Responses with Varying Gain (Motor and Load Inertia)')
% xlabel('Time(s)')
% ylabel('Angular Velocity (rad/s)')
% legend('Kp=1','Kp=2','Kp=5')
% figure(6)
% plot(t_msl1,R_msl1,'r-',t_msl2,R_msl2,'b--',t_msl5,R_msl5,'g-.')
% title('Step Responses with Varying Gain (Motor, Spring and Load Inertia)')
% xlabel('Time(s)')
% ylabel('Angular Velocity (rad/s)')
% legend('Kp=1','Kp=2','Kp=5')
% figure(7)
% plot(t_m2,R_m2,'r-',t_ml2,R_ml2,'b--',t_msl2,R_msl2,'g-.')
% title('Step Responses with Three Plant Configurations (Kp = 2)')
% xlabel('Time(s)')
% ylabel('Angular Velocity (rad/s)')
% legend('Motor alone','Motor and Load','Motor, Spring and Load')
