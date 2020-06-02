function [ ]=fourier_fft( )
%% 使用fft对任意激励进行傅里叶变换
wave=textread('.\GM1-TH.txt', '' ,'headerlines',1); % 读取地震信息
time=wave(:,1);
acc=wave(:,2);
dt=time(2)-time(1); % 采样间隔

% Tg=2;
% dt=0.001; % 采样间隔
% time=0:dt:20; % 激励的采样时刻
% acc=ag(Tg,time); % 变换前的激励

n=length(acc); % 样本点个数

df=1/time(end); % 频率间隔
f=0:df:(n-1)*df'; % 频率序列

%% FFT变换
y=fft(acc,n); % 作傅里叶变换
amp=abs(y); % 振幅（非真实振幅）
amp=amp/(n/2); % 真实振幅

%% 逆变换
sum=0;
for item=1:n
    sum=sum+y(item)*exp(i*2*pi*f(item)*time);
end

%% 绘图
close
subplot(2,1,1)
plot(f(1:n/2),amp(1:n/2))

subplot(2,1,2)
plot(time,acc)
hold on
plot(time,sum/n)

function [y]=ag(Tg,x)
%% 定义激励函数
Rem=rem(x,Tg);
Logic1=abs(Rem)>=0 & abs(Rem)<=Tg/2; % 判断大于等于零的位置
Logic2=~Logic1; % 判断小于零的位置
y=Logic1*1+Logic2*(-1);
