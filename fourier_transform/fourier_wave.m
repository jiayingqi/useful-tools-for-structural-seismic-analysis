function [ ]=fourier_wave()
%% 本函数用于求地震波的功率谱
clear
clc

%% 读取地震波
wave=textread('.\地震波\RSN12_KERN.PEL_PEL090.AT2', '' ,'headerlines',4); % 读取地震波

dt=0.005; % 地震波时间间隔(s)

wave=wave'; %行转列
wave=wave(:); %归为一列
wave=wave'*9.8; % m/s^2

%% 傅里叶变换
Omg=0:0.1:20; % 频域范围(Hz)
item=length(Omg)-1; % 谐波系数的项数
for n=-item:item % 求谐波系数
    item=item+1;
    cn(item)=1/Tg*trapz(t_range,wave.*exp(-i*n*Omg*t_range));
end

