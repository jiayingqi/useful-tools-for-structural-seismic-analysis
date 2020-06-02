function []=WaveJudgement()
%% 本函数用于判别单条地震波是否符合抗规要求

% 作者：贾英琦 王超
% 时间：2018/12/06
% 依据：《GB 50011-2010建筑抗震设计规范》(2016年版)
%      《星海诺富特高层框架-剪力墙酒店设计（大连理工大学本科毕业设计）》-贾英琦
% 说明：抗规对地震波的要求包括5个方面：（1）周期点误差；（2）底部剪力；（3）持时；
%     （4）数量；（5）峰值。本程序仅判别（1）周期点误差和（3）持时要求。

clear
clc
close all

%% 结构信息
T=[5.6 5.2 4.1]; % 高塔周期(s)
% T=[0.746 0.734 0.553]; % 矮塔周期(s)
ksi=0.05; % 阻尼比

%% 场地信息
switch 2
    case 1 % 根据场地信息进行计算
        group=1; % 场地分组： 1,2,3
        category=4; % 场地类别： 10,11,2,3,4
        impact='多遇'; % 设防烈度： '多遇','罕遇'
        intensity=7; % 设防烈度：6,7,7.5,8,8.5,9
        Tg=site(group,category); % 特征周期(s)
        alpha_max=situation(impact,intensity); % 影响系数峰值(g)
    case 2 % 直接赋值
        Tg=0.45;
        alpha_max=2.25; % 平台放大系数
end

%% 地震波信息
data_wave=textread('.\EQSignal-AW2-Acc.txt', '' , 'headerlines',0); % 读取地震波历程
data_spectra=textread('.\EQSignal-AW2-Acc-SP.txt', '' , 'headerlines',1); % 读取地震波反应谱

% time=data_wave(:,1);
dt=0.02;
time=(1:length(data_wave))*dt; % 时间序列

acc=data_wave(:,1); % 加速度历程
period=data_spectra(:,1); % 周期
alpha=data_spectra(:,2); % 影响系数

%% 周期点误差要求
T_range=0:0.1:6; % 周期序列
alphas_range=arrayfun(@(Ts)StandardSpectrum(Ts,ksi,alpha_max,Tg),T_range); % 标准反应谱

alphas_standard=arrayfun(@(Ts)StandardSpectrum(Ts,ksi,alpha_max,Tg),T) % 三个周期对应的标准曲线影响系数
alphas_calculate=arrayfun(@(Ts)CalculateAlpha(Ts,period,alpha),T) % 三个周期对应的地震波曲线影响系数

err=(alphas_calculate-alphas_standard)./alphas_standard
if sum(abs(err)<=repmat(0.20,[1,3]))==3
    disp('【满足】周期点误差要求')
else
    disp('【不满足】周期点误差要求')
end

%% 持时要求
item1=find(abs(acc)>0.1*max(abs(acc)),1,'first'); % 首次达到PGA的10%
item2=find(abs(acc)>0.1*max(abs(acc)),1,'last'); % 最后一次达到PGA的10%
tc=time(item2)-time(item1);
if tc>=max(15,5*max(T))
    disp('【满足】持时要求')
else
    disp('【不满足】持时要求')
end

%% 绘图
blue=[96 157 202]/256;
orange=[255 160 65]/256;
green=[56 194 93]/256;
pink=[255 91 78]/256;
purple=[184 135 195]/256;
gray=[164 160 155]/256;
fontsize=20;

figure('position',[100 100 1000 600])
subplot(2,1,1)
plot(T_range,alphas_range,'linewidth',3)
hold on
plot(period,alpha,'linewidth',3)
limit1=max(max(alphas_range),max(alpha));
plot([T(1),T(1)],[0,limit1*1.2],'linewidth',3)
plot([T(2),T(2)],[0,limit1*1.2],'linewidth',3)
plot([T(3),T(3)],[0,limit1*1.2],'linewidth',3)
set(legend('\alpha standard','\alpha calculate','\itT\rm1','\itT\rm2','\itT\rm3'),...
    'Fontname', 'Times New Roman','FontSize',fontsize,'EdgeColor',gray,'linewidth',1.5)
set(gca,'linewidth',2,'Fontname','Times New Roman','FontSize',fontsize)
set(title('Acceleration response spectrum'),'Fontname', 'Times New Roman','FontSize',fontsize)
set(xlabel('Period (s)'),'Fontname', 'Times New Roman','FontSize',fontsize)
set(ylabel('\alpha (g)'),'Fontname', 'Times New Roman','FontSize',fontsize)

subplot(2,1,2)
plot(time,acc,'linewidth',1.5)
hold on
limit2=max(abs(acc))*1.1;
plot([time(item1),time(item1)],[-limit2,limit2],'linewidth',3,'color',pink)
plot([time(item2),time(item2)],[-limit2,limit2],'linewidth',3,'color',pink)
plot([0,max(time)],[0.1*max(abs(acc)),0.1*max(abs(acc))],'linewidth',3,'color',pink)
plot([0,max(time)],[-0.1*max(abs(acc)),-0.1*max(abs(acc))],'linewidth',3,'color',pink)
set(gca,'linewidth',2,'Fontname','Times New Roman','FontSize',fontsize)
set(title('Acceleration history'),'Fontname', 'Times New Roman','FontSize',fontsize)
set(xlabel('Time (s)'),'Fontname', 'Times New Roman','FontSize',fontsize)
set(ylabel('Acc. (g)'),'Fontname', 'Times New Roman','FontSize',fontsize)

function alphas=CalculateAlpha(Ts,period,alpha)
alpha1=alpha(find(period>Ts,1,'first')); % 首次大于Ts时的alpha值
alpha2=alpha(find(period>Ts,1,'first')-1); % 首次小于Ts时的alpha值
alphas=(alpha1+alpha2)/2; % 取平均值

function alphas=StandardSpectrum(Ts,ksi,alpha_max,Tg)
% 计算标准反应谱
gamma=0.9+(0.05-ksi)/(0.3+6*ksi); % 曲线下降段的衰减指数
yta1=0.02+(0.05-ksi)/(4+32*ksi); % 直线下降段的下降斜率调整系数
yta2=1+(0.05-ksi)/(0.08+1.6*ksi); % 阻尼调整系数
if yta1<0
    yta1=0;
end
if yta2<0.55
    yta2=0.55;
end
if Ts<=0.1
    alphas=alpha_max*(0.45*(0.1-Ts)/0.1+yta2*Ts/0.1);
elseif Ts<=Tg
    alphas=yta2*alpha_max;
elseif Ts<=5*Tg
    alphas=(Tg/Ts)^gamma*yta2*alpha_max;
else
    alphas=(yta2*0.2^gamma-yta1*(Ts-5*Tg))*alpha_max;
end
    
function Tg=site(group,category)
% 确定场地特征周期
switch category
    case 10
        category=1;
    case 11
        category=2;
    case 2
        category=3;
    case 3
        category=4;
    case 4
        category=5;
end
data_site=...
    [0.20 0.25 0.35 0.45 0.65;
    0.25 0.30 0.40 0.55 0.75;
    0.30 0.35 0.45 0.65 0.90];
Tg=data_site(group,category);

function alpha_max=situation(impact,intensity)
% 确定影响系数峰值
switch impact
    case '多遇'
        impact=1;
    case '罕遇'
        impact=2;
end
switch intensity
    case 6
        intensity=1;
    case 7
        intensity=2;
    case 7.5
        intensity=3;
    case 8
        intensity=4;
    case 8.5
        intensity=5;
    case 9
        intensity=6;
end
data_alpha=...
    [0.04 0.08 0.12 0.16 0.24 0.32;
    0.28 0.50 0.72 0.90 1.20 1.40];
alpha_max=data_alpha(impact,intensity);