function [ ]=PowerSpectralDensityWelch()
%% 本函数使用韦尔奇功率谱密度估计法（加窗交叠法）求解地震波的功率谱
% 参考链接
% https://www.ilovematlab.cn/thread-492083-1-1.html
% https://wenku.baidu.com/view/5bcd1214fad6195f312ba641.html
clc
clear

%% 读取地震波信号
file=textread('.\GM6-TH.txt', '' , 'headerlines',1);
g=9.8; % 重力加速度 m/s^2
ug=file(:,2)*g;
dt=0.005; % 采样间隔 s
fs=1/dt; % 采样频率 Hz
N=length(ug); % 信号长度
t=(1:N)*dt; % 时间矢量


%% 求功率谱密度PSD
ug = ug - mean(ug); % 去除直流成分
win = hanning(512, 'periodic'); % 海明窗函数，返回一个512点的对称序列

[PSD, f] = pwelch(ug, win, 256, N, fs, 'twosided'); % 韦尔奇功率谱密度估计
% [PSD,F] = PWELCH(X,WINDOW,NOVERLAP,NFFT,Fs)
% PSD--------------------输入信号的功率谱密度的估计值
% f----------------------频率序列，当fs的单位为样本/s时，f的单位为Hz
% X----------------------输入信号
% WINDOW-----------------用于划分信号的窗口，信号会被划分为若干段，段数等于窗口长度
% NOVERLAP---------------重叠样本的长度，默认为窗口长度的50%
% N----------------------离散傅里叶变换DFT中点的数量
% fs---------------------采样频率
% 'twosided''onesided'---双边功率谱or单边功率谱，后者幅值是前者两倍

%% 方差
disp('频域求解的方差')
sum(PSD(1:end-1).*diff(f))

disp('时域求解的方差')
rms(ug)^2 % 均方根的平方
sum(ug.^2*dt)/dt/N; % 先对平方积分，再除以持时

%% 绘图
close all


figure('position',[100 100 900 700])
blue=[96 157 202]/256;
orange=[255 160 65]/256;
green=[56 194 93]/256;
pink=[255 91 78]/256;
purple=[184 135 195]/256;
gray=[164 160 155]/256;

fontsize=15;
subplot(2,1,1)
plot(t,ug,'linewidth',2,'color',orange)
title('Time history of acceleration','Fontname', 'Times New Roman','FontSize',fontsize)
set(xlabel('Time [s]'),'Fontname', 'Times New Roman','FontSize',fontsize)
set(ylabel('Acceleration [m/s^2]'),'Fontname', 'Times New Roman','FontSize',fontsize)
grid on

subplot(2,1,2)
plot(f,PSD,'linewidth',2,'color',green);
title('Power spectral density (PSD)','Fontname', 'Times New Roman','FontSize',fontsize)
set(xlabel('Frequency [Hz]'),'Fontname', 'Times New Roman','FontSize',fontsize)
set(ylabel({'Power spectral density (PSD)';'[m^2/s^3] or [(m/s^2)^2/Hz]'}),'Fontname', 'Times New Roman','FontSize',fontsize)
% xlim([0 10])
grid on