function [ ]=Wave2PSDF( ) 
%% 本函数用于求解地震加速度的功率谱密度函数
clear
clc

%% 地震信息输入
ug=textread('.\地震波\RSN12_KERN.PEL_PEL090.AT2', '' ,'headerlines',4); % 读取地震波
% ug=textread('.\地震波\RSN17_SCALIF_SLO234.AT2', '' ,'headerlines',4); % 读取地震波
% ug=textread('.\地震波\RSN1599_DUZCE_ATS030.txt', '' ,'headerlines',4); % 读取地震波

ug=ug'; %行转列
ug=ug(:); %归为一列
dt=0.02; % 地震波采样间隔（s）
f=1/dt; % 采样频率（Hz）
g=9.8; % 重力加速度（m/s^2）
ug=ug*g; % 真实地震历程

ng=length(ug); % 当前样本点数目
time=linspace(dt,ng*dt,ng);

n=nextpow2(ng); % 本函数用于求解满足方程2^n≥|ng|的n的最小值
m=2^n; % 此n满足n≥|ng|且是以2为底数的幂函数，为扩展后的样本点的数目
% 当样本容量满足2的n次幂时，才可使用快速傅里叶变换

%% 计算单边功率谱密度函数
Y=fft(ug,m);       
% 表示利用快速傅里叶变换的方法计算ug向量的离散傅里叶变换（DFT），变换后的维数为n，不足则取0，冗余则截断
Sg=Y.*conj(Y)/m; % 功率谱密度（m^2/s^3）
P=Y.*conj(Y); % 功率谱（m^2/s^3）
Amp=abs(Y); % 傅里叶谱（m/s^2）
% 依据：钟骁珺, 基于MATLAB的地震波动力特性分析, 2010.

f_range=(1:m/2)/m*f; % 频率序列（Hz）
f_range=f_range';
Sg=Sg(1:m/2);
P=P(1:m/2);
Amp=Amp(1:m/2);

%% 绘图
close all
blue=[96 157 202]/256;
orange=[255 160 65]/256;
green=[56 194 93]/256;
pink=[255 91 78]/256;
purple=[184 135 195]/256;
gray=[164 160 155]/256;
fontSize=12;
item=500; % 绘制项数

subplot(2,2,1)
plot(f_range(1:item),Sg(1:item),'linewidth',1.5,'color',blue)
set(xlabel('Frequency \itf \rm(Hz)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(ylabel('PSDF of Acc. (m^2/s^3)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(gca,'Fontname', 'Times New Roman','FontSize',fontSize)

subplot(2,2,2)
plot(time,ug,'linewidth',1.5,'color',blue)
set(xlabel('Time \itt \rm(s)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(ylabel('Acc. (m/s^2)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(gca,'Fontname', 'Times New Roman','FontSize',fontSize)

subplot(2,2,3)
plot(f_range(1:item),P(1:item),'linewidth',1.5,'color',orange)
set(xlabel('Frequency \itf \rm(Hz)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(ylabel('Power spectrum of Acc. (m^2/s^3)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(gca,'Fontname', 'Times New Roman','FontSize',fontSize)

subplot(2,2,4)
plot(f_range(1:item),Amp(1:item),'linewidth',1.5,'color',green)
set(xlabel('Frequency \itf \rm(Hz)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(ylabel('Fourier spectrum of Acc. (m/s^2)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(gca,'Fontname', 'Times New Roman','FontSize',fontSize)

set(gcf,'position',[200,200,900,700])