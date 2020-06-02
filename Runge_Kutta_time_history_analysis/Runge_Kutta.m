function [u,du] = Runge_Kutta(wave,dt,n,M,C,K,E)
N=length(M); % 自由度数量
% 对位移、速度、加速度和地面加速度增量矩阵赋初值0
u(:,1)=zeros(N,1); 
du(:,1)=zeros(N,1);
ddu(:,1)=-E*wave(1);

P=M^-1*C; % M^-1表示对矩阵M求逆
Q=M^-1*K;

for i=1:(n-1)
    % 四阶Runge_Kutta参数
    L0=dt*(-P*du(i)-Q*u(i)-wave(i));
    L1=dt*(-P*(du(i)+L0/2)-Q*(u(i)+dt*du(i)/2)-(wave(i)+wave(i+1)/2));
    L2=dt*(-P*(du(i)+L1/2)-Q*(u(i)+dt*du(i)/2++dt*L0/4)-(wave(i)+wave(i+1)/2));
    L3=dt*(-P*(du(i)+L2)-Q*(u(i)+dt*du(i)++dt*L1/2)-wave(i+1));
    
    % 位移增量
    delta_u=dt*du(i)+dt*(L0+L1+L2)/6;
    delta_du=(L0+2*L1+2*L2+L3)/6;
    %求出迭代后各个变量的值
    u(:,i+1)=u(:,i)+delta_u;
    du(:,i+1)=du(:,i)+delta_du;  
end