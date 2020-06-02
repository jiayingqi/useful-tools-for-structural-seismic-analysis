%单自由度主结构时程分析
function [d,v,a] = SDOF_time_history_Newmark(m,c,k,acc,dt)

    long=length(acc);         %用于计算输入地震波的长度
    time=0:dt:(long-1)*dt;
    
    %由于TVMD自身有质量，等同于自由度多1倍
    F=-m*acc;
    
    %定义初始位移，速度和加速度
    y(:,1)=zeros(1,1);      %初始位移
   dy(:,1)=zeros(1,1);      %初始速度
  ddy(:,1)=zeros(1,1);      %初始加速度
 
  beta=1/4;
 for n=2:long;
    
    %计算地面运动加速度增量
	zacc=acc(n)-acc(n-1);
    
    dF_=-m.*zacc+m*(1./beta.*dy(:,n-1)./dt+1/2./beta.*ddy(:,n-1))+1/2./beta*c*dy(:,n-1);
    dk_=k+2.*c./dt+4.*m./(dt.^2);
    
     %计算响应位移增量
    zy(:,n)=dk_\dF_;        
    %计算响应位移
    y(:,n)=y(:,n-1)+zy(:,n); 
    
    %计算响应速度增量
    zdy(:,n)=1/2./beta.*zy(:,n)./dt-1/2./beta.*dy(:,n-1);
    %计算响应速度
    dy(:,n)=dy(:,n-1)+zdy(:,n);
 
    %加速度响应
    Fs(:,n)=k*y(:,n);
    ddy(:,n)=m\(F(:,n)-c*dy(:,n)-Fs(:,n));
    
	
 end
a=ddy;            %相对加速度
v=dy;             %相对速度
d=y;              %相对位移
