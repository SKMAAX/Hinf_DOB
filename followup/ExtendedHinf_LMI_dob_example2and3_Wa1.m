%% 拡張H∞制御のLMI解法で外乱オブザーバ設計
% 「H∞制御と外乱オブザーバ」例題2&例題3(Wa1=0.1)への適用
%  結果が微妙に合わない

clear;
s = tf('s');
%% 
% 相補感度関数Taに対する重みWaの設定

Da = 0.1;
Wa = tf([Da],[1]);
% bodemag(Wa);
%% 
% 制御対象の記述

A0 = [-1];
b = [1]; 
c = [1]; 
[n0,~] = size(A0);

%% 
% 一般化制御対象の記述

% A = [A0 b zeros(n0,na);zeros(1,n0) 0 zeros(1,na); zeros(na,n0) zeros(na,1) Aa];
A = [A0 b; 0 0];
[n,~] = size(A);
B1 = [zeros(n0,2);1 0];
B2 = [1;0];
C1 = [zeros(1,n0) 1];
C2 = [c 0];
D12 = [1];
D21 = [0 Da];
D12_dag = 1;
D21_dag = [0;1/Da];

E12 = D12'*D12;
E21 = D21*D21';

P = ss(A,B2,C2,0);
Ps = tf(P(1));
bode(Ps);
%% 
% LMIで使用する各行列の定義

U0 = zeros(n,n);
U1 = eye(n);
U1_dag = pinv(U1);
A_U = U1*A*U1_dag;
C1_U = C1*U1_dag;
C2_U = C2*U1_dag;
B1_U = U1*B1;
U_L = [];
%% 
% オブザーバゲインの算出

% SDPソルバとLMIパーサのパスの設定
addpath(genpath('～～sedumi-masterへのパス～～'))
addpath(genpath('～～YALMIP-masterへのパス～～'))

Y_U=sdpvar(n,n,'sy');
% X_V=sdpvar(n-1,n-1,'sy');
X_V = 0;
M=sdpvar(1,n-1); % M=F*V_1*X_V;
N=sdpvar(n,1); % N=Y_U*U_1*L;
gamma = sdpvar(1);

% -------------------
LMI=[]; % LMI initialized

ep = 1e-6;
LMI=[LMI, Y_U>=ep*eye(n)];
LMI2 = [A_U'*Y_U+Y_U*A_U+C2_U'*N'+N*C2_U+C1_U'*C1_U   Y_U*B1_U+N*D21;
    B1_U'*Y_U+D21'*N'  -eye(n)*gamma];
LMI=[LMI, LMI2<=-ep*eye(2*n)];

% -------------------
% ops=sdpsettings; ops.shift=1e-6;
sol=solvesdp(LMI,gamma);
% sol=solvesdp(LMI);

if sol.problem~=1;
	Y_U_a = double(Y_U)
	N_a = double(N)
    gamma_a = double(gamma)
	% -------------------
	format short e
	pres = checkset(LMI)
else
    disp(sol.info)
end
Y = U1_dag * inv(Y_U_a) * U1_dag';

L = inv(U1) * inv(Y_U_a)*N_a
F = -D12_dag*C1
%% 
% 制御器K(s)の算出

K = ss(A+B2*F+L*C2,[L B2],[-F;C2],[0 -1;1 0]);
tfK = tf(K);
Ks = tfK(1,1);
%% 
% 感度関数と相補感度関数

S = tf(feedback(1,Ks*Ps,+1)); %inv(1-P*K);
T = tf(feedback(Ks*Ps,1,+1)); %inv(1-P*K)*P*K;
Ta = tf(feedback(Ks,Ps,+1)); %inv(1-P*K)*K;
%% 
% 
% 
% 相補感度関数Tと、重み関数Waの大小関係確認

bodemag(Ta,1/Wa)
title('準相補感度関数Taと、重みWaの大小関係確認')
legend('Ta','1/Wa')

l1 = L(1)
l2 = L(2)
Stemp = tf([1 1+l1],[1 l1+1 l2]);
Ttemp = tf([l2 l2],[1 1+l1 l2]);
bodemag(Ta,Ttemp)