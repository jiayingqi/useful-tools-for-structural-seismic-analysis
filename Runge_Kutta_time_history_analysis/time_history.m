function []=time_history()
%% 用Newmark法对结构体系进行时程分析
clear
clc
%% 输入结构信息
m=2000; %原结构质量，kg
omega=2*pi/1.5; %原结构频率，rad/s
k=m*omega^2; %原结构刚度，N/m
ksi=0.03; %原结构阻尼比
c=2*ksi*omega*m; %原结构阻尼系数，N·s/m

%% 矩阵求解
[M,C,K,E] = matrix_shear_building(m, c, k);

%% 读取激励
wave=textread('.\地震波\RSN12_KERN.PEL_PEL090.AT2', '' ,'headerlines',4); % 读取地震波

dt=0.005; % 地震波时间间隔(s)

wave=wave'; %行转列
wave=wave(:); %归为一列
wave=wave'*9.8; % m/s^2
n=length(wave);

%% Newmark求解位移时程响应
[u,du,ddu] = Newmark_belta(wave,dt,n,M,C,K,E);
[u1,du1] = Runge_Kutta(wave,dt,n,M,C,K,E);

t=linspace(0.005,n*0.005,n);

%% 后处理
close all
green=[64/256,116/256,52/256];
blue=[7/256,151/256,237/256];
orange=[248/256,147/256,29/256];
red=[219/256,69/256,32/256];
brown=[107/256,90/256,92/256];
gray=[151/256,151/256,151/256];

figure('position',[300,300,1200,600]) % 位移

subplot(2,1,1)
plot(t,u,'linewidth',2,'color',orange)
set(xlabel('Time (s)'),'Fontname', 'Times New Roman','FontSize',15)
set(ylabel('Displacement (m)'),'Fontname', 'Times New Roman','FontSize',15)
set(legend('Newmark-\beta'),'Fontname', 'Times New Roman','FontSize',15)
set(gca,'Fontname', 'Times New Roman','FontSize',15)

subplot(2,1,2)
plot(t,u1,'linewidth',2,'color',blue)
set(xlabel('Time (s)'),'Fontname', 'Times New Roman','FontSize',15)
set(ylabel('Displacement (m)'),'Fontname', 'Times New Roman','FontSize',15)
set(legend('Runge-Kutta'),'Fontname', 'Times New Roman','FontSize',15)
set(gca,'Fontname', 'Times New Roman','FontSize',15)
