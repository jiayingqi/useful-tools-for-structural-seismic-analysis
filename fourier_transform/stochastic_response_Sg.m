function  [Omega, Sx, Sigma_X, Sigma_XP]=stochastic_response_Sg(lamda1, Phi1, r, Sg, Omega)

Omega=Omega';
N=length(lamda1)/2; % 自由度个数（包含惯容自由度）
n=length(Omega); % 功率谱密度的长度

% Omega=Omega';

%% 求Z向量
Z=ones(n,1)*r';
Z=Z'./(Omega*sqrt(-1)-lamda1);
Sg=Sg';
Y=Phi1*Z.*sqrt(Sg);

%% 求自谱
X=Y(N+1:2*N,:); % 取Y的下半部分（位移传递函数）
XP=Y(1:N,:); % 取Y的上半部分（速度传递函数）
Sx=conj(X).*X; % 位移自谱
Sxp=conj(XP).*XP; % 速度自谱
Sx=1*Sx; % 单边功率谱密度
Sxp=1*Sxp;

%% 求位移和速度
dOmega=Omega(2:n)-Omega(1:n-1);
Sum1=dOmega.*Sx(:,2:n);
Sigma_X=sqrt(sum(Sum1,2)); %sum(x,2)表示对x矩阵的每一行分别求和，并存入列向量
Sum2=dOmega.*Sxp(:,2:n);
Sigma_XP=sqrt(sum(Sum2,2));
