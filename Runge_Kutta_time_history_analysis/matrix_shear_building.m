function [M,C,K,E] = matrix_shear_building(m,c,k)
zeta=0.05;
omg = sqrt(k./m);
N=length(m);
M=diag(m);

%% ×èÄá¾ØÕó
for i=1:N-1
    Cd(i,i)=c(i)+c(i+1);
    Cd(i,i+1)=-c(i+1);
    Cd(i+1,i)=Cd(i,i+1);
end
Cd(N,N)=c(N);

%% ¸Õ¶È¾ØÕó
for i=1:N-1
    K(i,i)=k(i)+k(i+1);
    K(i,i+1)=-k(i+1);
    K(i+1,i)=K(i,i+1);
end
K(N,N)=k(N);

%% ÇóÆµÂÊºÍ×èÄá
% omega2 = eig(K, M);
% omega = sqrt(omega2);
% 
% w_i = omega(1);
% w_j = omega((floor(0.618*N))); % Ä£Ì¬½Ø¶Ï
% 
% a0 = zeta*2.0*w_i*w_j/(w_i+w_j);
% a1 = zeta*2.0/(w_i+w_j);

% Cs = a0*M + a1*K % ÈğÀû×èÄá
Cs = 0;
C = Cs + Cd; % È«²¿×èÄá
E=ones(N,1);
