function [ ] = Parseval( )
%% 本程序用于验证Parseval’s theorem（时域和频域的能量相等）
% 该定理的数学表达形式为sigma^2=integral(0,Tg,x(t)^2/Td)=integral(-∞,+∞,S(omega))
% sigma^2为函数x(t)的方差，x(t)应满足均值为0的条件，S(omega)为x(t)的功率谱密度函数，Tg为x(t)的持时
% 参考：https://blog.csdn.net/enjoy_learn/article/details/60978606

clear
clc

%% 地震动信息输入
switch 2
    case 1 % 人为定义
        dt = 0.001; % 采样间隔，s
        fs = 1/dt; % 采样频率，Hz
        Td = 30;
        t = 0:dt:Td; % 时间序列
        N = length(t); % 采样点数
        fx = 15*cos(0.5*pi*t)+8*cos(1*pi*t); % 样本函数
    case 2 % 读取已有地震动
        dt = 0.005;
        fs = 1/dt; % 采样频率，Hz
%         fx=textread('.\地震波\RSN17_SCALIF_SLO234.AT2', '' ,'headerlines',4); % 读取地震波
%         fx=textread('.\地震波\RSN1599_DUZCE_ATS030.txt', '' ,'headerlines',4); % 读取地震波
%         fx=textread('.\地震波\RSN1540_CHICHI_TCU115-N.AT2', '' ,'headerlines',4); % 读取地震波
%         fx=textread('.\地震波\RSN718_SUPER.A_A-IVW090.AT2', '' ,'headerlines',4); % 读取地震波
        fx=textread('.\地震波\RSN16_NCALIF.AG_A-FRN044.AT2', '' ,'headerlines',4); % 读取地震波

        fx = fx'; % 行转列
        g = 9.8; % m/s^2
        fx = fx(:)*g;
        N = length(fx); % 采样点数
        Td = (N-1)*dt;
        t = 0:dt:Td; % 时间序列
end

%% 时域求解激励方差
disp('时域角度求激励方差')
variance_time = sum(fx.^2)/Td*dt % 时域角度求方差

%% 频域求解激励方差
n=nextpow2(N); % 本函数用于求解满足方程2^n≥|N|的n的最小值
m=2^n; % 扩展样本点数目，此n满足n≥|ng|且是以2为底数的幂函数，为扩展后的样本点的数目
% 当样本容量满足2的n次幂时，才可使用快速傅里叶变换

window = hanning(N); % 汉宁窗函数
[Pxx_period,f_period] = periodogram(fx,[],m,fs); % periodogram和pwelch可以用于PSD的估计，如有必要，可以将第二项改为window
df = f_period(2)-f_period(1); % 输出频率的间隔，不同于采样频率fs
S = Pxx_period/m/dt/df; % 考虑幅值修正的PSD

disp('频域角度求激励方差')
variance_frequency = sum(S)*df % 频域角度求方差

%% 时域角度求解响应方差
m = 200e3; % 原结构质量，kg
Ts = 1.5; % 原结构周期，s
omega = 2*pi/Ts; % 原结构频率，rad/s
k = m*omega^2; % 原结构刚度，N/m
zeta = 0.03; % 原结构阻尼比
c = 2*zeta*omega*m; %原结构阻尼系数，N·s/m

[u,du,ddu] = Newmark_belta(fx,dt,N,m,c,k,1);
disp('时域角度求响应方差')
variance_time_response = sum(u.^2)/Td*dt
% variance_time_response = rms(u)^2

%% 频域角度求响应方差
[lamda, Phi, r] = complex_modes(m,c,k,1);
S = S/(2*pi); % 考虑幅值修正的PSD，对的频率单位为rad/s
Omega = f_period*2*pi; % 圆频率，单位为rad/s
[Omega, Sx, Sigma_X, Sigma_XP] = stochastic_response_Sg(lamda, Phi, r, S, Omega);
Sx = Sx*2*pi; % 将S(omega)转未S(f)
disp('频域角度求响应方差')
variance_frequency_response = Sigma_X^2

%% 绘图
close all
blue=[96 157 202]/256;
orange=[255 160 65]/256;
green=[56 194 93]/256;
pink=[255 91 78]/256;
purple=[184 135 195]/256;
gray=[164 160 155]/256;
fontSize=15;

set(gcf,'position',[200,200,900,700])
subplot(2,2,1)
plot(t,fx,'-','linewidth',2,'color',blue)
set(xlabel('Time \itt \rm(s)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(ylabel('Acc. \ita \rm(m/s^2)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(gca,'Fontname', 'Times New Roman','FontSize',fontSize)

subplot(2,2,2)
plot(f_period,S,'-','linewidth',3,'color',pink)
set(xlabel('Frequency \itf \rm(Hz)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(ylabel('PSDF of Acc. \itS \rm(m^2/s^3)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(gca,'Fontname', 'Times New Roman','FontSize',fontSize)
xlim([0,20])

subplot(2,2,3)
plot(t,u,'-','linewidth',2,'color',green)
set(xlabel('Time \itt \rm(s)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(ylabel('Disp. response \itu \rm(m)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(gca,'Fontname', 'Times New Roman','FontSize',fontSize)

subplot(2,2,4)
plot(f_period,Sx,'-','linewidth',3,'color',pink)
set(xlabel('Frequency \itf \rm(Hz)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(ylabel('PSDF of Acc. \itS \rm(m^2/s^3)'),'Fontname', 'Times New Roman','FontSize',fontSize)
set(gca,'Fontname', 'Times New Roman','FontSize',fontSize)
xlim([0,20])