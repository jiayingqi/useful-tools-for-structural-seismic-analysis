function [ ]=test()
clear
clc

%% 本函数用于验证傅里叶变换
m=2000; %原结构质量，kg
omega=2*pi/1.5; %原结构频率，rad/s
k=m*omega^2; %原结构刚度，N/m
Tg=4; % 周期激励的周期
1/k
Omg=2*pi/Tg; % 周期激励的频率
dt=0.01; % 采样间隔
t=0:dt:20; % 激励的采样时刻
y=ag(Tg,t);
item=1; % 傅里叶变换的项数

t_range=0:0.01:Tg; % 求谐波系数时的积分序列
p0=1/Tg*trapz(t_range,ag(Tg,t_range));
for i=1:item % 求谐波系数
    pc(i)=2/Tg*trapz(t_range,ag(Tg,t_range).*cos(i*Omg*t_range));
    ps(i)=2/Tg*trapz(t_range,ag(Tg,t_range).*sin(i*Omg*t_range));
end
digits(4) % 精度控制
rn=(1:item)*Omg/sqrt(k/m);

%% 求变换后的激励和响应
P=p0; % 激励
U=p0/k; % 响应
% for i=1:item
%     P=P+pc(i)*cos(i*Omg*t)+ps(i)*sin(i*Omg*t);
%     U=U+(pc(i)*cos(i*Omg*t)+ps(i)*sin(i*Omg*t))/k/(1-rn(i)^2);
% end
% U1=1/k/(1-Omg^2/omega^2)*(sin(Omg*t)-Omg/omega*sin(omega*t));
% U2=1/k/(1-Omg^2/omega^2)*(cos(Omg*t)-cos(omega*t));
U3=1/k*(1-cos(omega*t));

%% 求变换前的响应
[u,du,ddu] = Newmark_belta(y,dt,length(y),m,0,k,1);

%% 画图验证
close all
figure('position',[100 100 700 400])
subplot(4,1,1)
plot(y)
hold on
plot(P)
ylim([-2,2])
legend('变换前','变换后')

subplot(4,1,2)
% plot(u)
% hold on
plot(U*m)
% legend('变换前','变换后')


subplot(4,1,3)
plot(U3*m)

subplot(4,1,4)
plot(-u)

function [y]=ag(Tg,x)
%% 定义激励函数
% y=cos(2*pi/Tg*x);
y=ones(1,length(x));