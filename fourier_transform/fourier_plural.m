function []=fourier_plural()
%% 本函数用于验证复数形式的傅里叶变换
% 本函数中考虑了结构阻尼（这在实数形式的傅里叶变换也可以考虑，但实数形式无法考虑任意激励）
clear
clc

%% 结构信息
m=2000; %原结构质量，kg
omega=2*pi/1.5; %原结构频率，rad/s
k=m*omega^2; %原结构刚度，N/m
ksi=0.03; %原结构阻尼比
c=2*ksi*omega;
Tg=4; % 周期激励的周期ega*m; %原结构阻尼系数，N·s/m

%% 激励信息
Omg=2*pi/Tg; % 周期激励的频率
dt=0.01; % 采样间隔
t=0:dt:50; % 激励的采样时刻
y=ag(Tg,t); % 加速度序列
F=-m*ag(Tg,t); % 荷载序列
item=5; % 傅里叶变换的项数

t_range=0:0.01:Tg; % 求谐波系数时的积分序列

item_cn=0; % 谐波系数的项数
for n=-item:item % 求谐波系数
    item_cn=item_cn+1;
    cn(item_cn)=1/Tg*trapz(t_range,-m*ag(Tg,t_range).*exp(-i*n*Omg*t_range)); % 前乘-m是将加速度转变为力
end

%% 变换后的荷载和响应
item_cn=0;
P=0;
U=0;
u_free=0;
omegad=omega*sqrt(1-ksi^2);
rn=(-item:item)*Omg/omega;
H=1./(1-rn.^2+2*ksi*rn*i)/k;
B=(-c-sqrt(c^2+4*m*k))/2/m;

for n=-item:item
    item_cn=item_cn+1;
    P=P+cn(item_cn)*exp(i*n*Omg*t); % 变换后的荷载
    U=U+cn(item_cn)*H(item_cn)*exp(i*n*Omg*t); % 变换后的响应
end


%% 变换前的响应
[u1,~,ddu] = Newmark_belta(y,dt,length(y),m,c,k,1); % 自己的程序
[u2,~,~] = SDOF_time_history_Newmark(m,c,k,y,dt); % 张老师的程序

%% 画图验证
close all
figure('position',[100 100 700 400])
subplot(3,1,1)
plot(t,F,'linewidth',2)
hold on
plot(t,P,'linewidth',2)
legend('变换前的地震荷载','变换后的地震荷载')

subplot(3,1,2)
plot(t,u1,'linewidth',2)
hold on
plot(t,U,':','linewidth',2)
legend('变换前的响应（时程分析）','变换后的响应（未考虑自由振动）')

subplot(3,1,3)
plot(t,U-u1,'linewidth',2)
legend('自由振动')

function [y]=ag(Tg,x)
%% 定义激励函数
Rem=rem(x,Tg);
Logic1=abs(Rem)>=0 & abs(Rem)<=Tg/2; % 判断大于等于零的位置
Logic2=~Logic1; % 判断小于零的位置
y=Logic1*1+Logic2*(-1);