function [u,du,ddu] = Newmark_belta(wave,dt,n,M,C,K,E)
N=length(M); % 自由度数量
% 对位移、速度、加速度和地面加速度增量矩阵赋初值0
u(:,1)=zeros(N,1); 
du(:,1)=zeros(N,1);
ddu(:,1)=-ones(N,1)*wave(1);
delta_ddy(:,1)=zeros(N,1);
for i=1:(n-1)
    % 地面运动加速度增量
    delta_ddy(:,1)=wave(i+1)-wave(i);
    % 等效刚度
    Ke=K+2*C/dt+4*M/dt^2;
    % 荷载增量
    delta_P=M*(-E.*delta_ddy+4*du(:,i)/dt+2*ddu(:,i))+2*C*du(:,i);      
    %位移增量
    delta_u=Ke^(-1)*delta_P;
    %速度增量
    delta_du=2*delta_u/dt-2*du(:,i);
    %加速度增量
    delta_ddu=4*delta_u/dt^2-4*du(:,i)/dt-2*ddu(:,i);
    %求出迭代后各个变量的值
    u(:,i+1)=u(:,i)+delta_u;
    du(:,i+1)=du(:,i)+delta_du;
    ddu(:,i+1)=ddu(:,i)+delta_ddu;
end