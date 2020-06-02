function [ ]=TestTimeFrequency()
%% 本函数用于验证时域和频域的吻合性
clear
clc

%% 结构信息
m=200; % 原结构质量，kg
T=1.5; % 原结构周期,s
omega=2*pi/T; % 原结构频率，rad/s
k=m*omega^2; % 原结构刚度，N/m
zeta=0.03; % 原结构阻尼比
c=2*zeta*omega*m; %原结构阻尼系数，N·s/m

%% 读取激励
ug=textread('F:\D-SMA-TMD\地震波\GM2-TH.txt', '' ,'headerlines',1); % 读取地震波
ug=ug(:,2);
ug=ug(:); % 归为一列
ug=ug'; 
ug=ug(1:1000);
dt=0.02;
n=length(ug);
t=linspace(dt,n*dt,n);
[Sg,Omega]=Wave2PSDF(ug,dt);

%% 频域响应求解
[M,C,K,E] = matrix_shear_building(m, c, k);
[lamda, Phi, r] = complex_modes(M,C,K,E);
[~, Sx, Sigma_X, Sigma_XP] = stochastic_response_Sg(lamda, Phi, r, Sg, Omega);

disp('频域求解的位移均方根')
Sigma_X

%% 时域响应求解
[u,du,ddu] = Newmark_belta(ug,dt,length(ug),m,c,k,1);
disp('时域求解的位移均方根')
rms(u)


close all
subplot(5,1,1)
plot(Omega,Sg)
subplot(5,1,2)
plot(t,ug)
subplot(5,1,3)
plot(t,ddu)
subplot(5,1,4)
plot(t,du)
subplot(5,1,5)
plot(t,u)
set(gcf,'position',[200,200,900,700])

function [Sg,Omega]=Wave2PSDF(ug,dt) 
%% 地震信息输入
% ug--------------------加速度序列
% dt--------------------地震波采样间隔（s）
ug=ug'; %行转列
ug=ug(:); %归为一列
f=1/dt; % 采样频率（Hz）

ng=length(ug); % 当前样本点数目
time=linspace(dt,ng*dt,ng);

n=nextpow2(ng); % 本函数用于求解满足方程2^n≥|ng|的n的最小值
m=2^n; % 此n满足n≥|ng|且是以2为底数的幂函数，为扩展后的样本点的数目
% 当样本容量满足2的n次幂时，才可使用快速傅里叶变换

%% 计算单边功率谱密度函数
Y=fft(ug,m);       
% 表示利用快速傅里叶变换的方法计算ug向量的离散傅里叶变换（DFT），变换后的维数为n，不足则取0，冗余则截断
Sg=Y.*conj(Y)/m; % 功率谱密度（m^2/s^3）
Sg=abs(Y)/m;
% 依据：钟骁珺, 基于MATLAB的地震波动力特性分析, 2010.
f_range=(1:m/2)/m*f; % 频率序列（Hz）
f_range=f_range';
Sg=Sg(1:m/2);
Omega=f_range*2*pi;