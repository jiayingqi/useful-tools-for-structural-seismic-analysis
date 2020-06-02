function [ ]=fourier_real_number( )
%% 本函数用于验证实数形式的傅里叶变换
clear
clc

%% 结构信息
m=2000; %原结构质量，kg
omega=2*pi/1.5; %原结构频率，rad/s
k=m*omega^2; %原结构刚度，N/m
Tg=4; % 周期激励的周期

%% 激励信息
Omg=2*pi/Tg; % 周期激励的频率
dt=0.01; % 采样间隔
t=0:dt:50; % 激励的采样时刻
y=ag(Tg,t); % 加速度序列
F=-m*ag(Tg,t); % 荷载序列
item=20; % 傅里叶变换的项数

%% 求谐波系数
t_range=0:0.01:Tg; % 求谐波系数时的积分序列
p0=1/Tg*trapz(t_range,-m*ag(Tg,t_range));
for i=1:item 
    pc(i)=2/Tg*trapz(t_range,-m*ag(Tg,t_range).*cos(i*Omg*t_range));
    ps(i)=2/Tg*trapz(t_range,-m*ag(Tg,t_range).*sin(i*Omg*t_range));
end
digits(4) % 精度控制
rn=(1:item)*Omg/sqrt(k/m);

%% 求变换后的激励和响应
P=p0; % 激励
U=p0/k; % 响应
U_free=p0/k*(1-cos(omega*t));
for i=1:item
    P=P+pc(i)*cos(i*Omg*t)+ps(i)*sin(i*Omg*t);
    U=U+1/k/(1-rn(i)^2)*(pc(i)*cos(i*Omg*t)+ps(i)*sin(i*Omg*t));
    U_free=U_free+1/k/(1-rn(i)^2)*(pc(i)*cos(i*Omg*t)+ps(i)*sin(i*Omg*t)+...
        pc(i)*cos(omega*t)-ps(i)*i*Omg/omega*sin(omega*t));
end

%% 求变换前的响应
[u,du,ddu] = Newmark_belta(y,dt,length(y),m,0,k,1);

%% 画图验证
close all
figure('position',[100 100 700 400])
subplot(2,1,1)
plot(t,F,'linewidth',2)
hold on
plot(t,P,'linewidth',2)
legend('变换前的地震荷载','变换后的地震荷载')

subplot(2,1,2)
plot(t,u,'linewidth',2)
hold on
plot(t,U,'-.','linewidth',2)
plot(t,U_free,':','linewidth',2)
legend('变换前的响应（时程分析）','变换后的响应（未考虑自由振动）','变换后的响应（考虑自由振动）')

function [y]=ag(Tg,x)
%% 定义激励函数
Rem=rem(x,Tg);
Logic1=abs(Rem)>=0 & abs(Rem)<=Tg/2; % 判断大于等于零的位置
Logic2=~Logic1; % 判断小于零的位置
y=Logic1*1+Logic2*(-1);

% y=sin(2*pi/Tg*x);

% y=2*ones(1,length(x));