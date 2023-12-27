%% 拡張H∞制御のLMI解法の数値例確認
% 

clear;
s = tf('s');

%% 
% 制御対象の記述

A = [0 1;0 -1];
[n,~] = size(A);
B1 = [0;1];
B2 = [0;1];
C1 = [1 0; 0 0];
C2 = [0 1];
D12 = [0; 0.1];
D21 = [0];

G = ss(A,[B1 B2],[C1;C2],[zeros(2,1) D12; D21 0]);
P = ss(A,B1,C2,0);
Ps = tf(P(1));
%% 
% LMIで使用する各行列の定義


V0 = zeros(2,2);
V1 = eye(2);
V1_dag = pinv(V1);
A_V = V1_dag*A*V1;
B1_V = V1_dag*B1;
B2_V = V1_dag*B2;
C1_V = C1*V1;
V_F = zeros(1,2);


U0 = [1 0];
U1 = [0 1];
U1_dag = pinv(U1);
A_U = U1*A*U1_dag;
C1_U = C1*U1_dag;
C2_U = C2*U1_dag;
B1_U = U1*B1;
U_L =-1;

%% 
% LMIによるゲインの算出

% SDPソルバとLMIパーサのパスの設定
addpath(genpath('～～sedumi-masterへのパス～～'))
addpath(genpath('～～YALMIP-masterへのパス～～'))

Y_U=sdpvar(1,1,'sy');
X_V=sdpvar(2,2,'sy');
M=sdpvar(1,2); % M=F*V_1 * X_V;
N=sdpvar(1,1); % N=Y_U*U_1*L;
% -------------------
ep = 1e-6;
LMI=[]; % LMI initialized

LMI=[LMI, Y_U>=ep*eye(length(Y_U))];
LMI=[LMI, X_V>=ep*eye(length(X_V))];

LMI1 = [A_V*X_V+X_V*A_V'+B2_V*M+M'*B2_V'+B1_V*B1_V'   X_V*C1_V'+M'*D12';
    C1_V*X_V+D12*M  -eye(2)];
LMI=[LMI, LMI1<=-ep*eye(length(LMI1))];

LMI2 = [A_U'*Y_U+Y_U*A_U+C2_U'*N'+N*C2_U+C1_U'*C1_U   Y_U*B1_U+N*D21;
    B1_U'*Y_U+D21'*N'  -eye(1)];
LMI=[LMI, LMI2<=-ep*eye(length(LMI2))];

LMI3 = [X_V V1_dag*U1_dag;
    U1_dag'*V1_dag' Y_U];
LMI=[LMI, LMI3>=-ep*eye(length(LMI3))];

% -------------------
% ops=sdpsettings; ops.shift=1e-6;
sol=solvesdp(LMI);

if sol.problem~=1;
	X_V_a = double(X_V)
	Y_U_a = double(Y_U)
	M_a = double(M)
	N_a = double(N)
	% -------------------
	format short e
	pres = checkset(LMI)
end
X = V1_dag' * inv(X_V_a) * V1_dag;
Y = U1_dag * inv(Y_U_a) * U1_dag';
F = [M_a*inv(X_V_a)] * inv([V1]) ;
L = inv([U0; U1]) * [U_L; inv(Y_U_a)*N_a];

%% 
% 制御器K(s)の算出

K = ss(A+B2*F+L*C2,[L B2],[-F;C2],[0 -1;1 0]);
tfK = tf(K);
Ks = tfK(1,1);
%% 
% 感度関数と相補感度関数

S = tf(feedback(1,Ps*Ks,+1)); %inv(1-P*K);
Ta = tf(feedback(Ks,Ps,+1)); %inv(1-P*K)*K;

pole(S) %閉ループ系の固有値
%% 
% w→z1、w→z2の閉ループ系伝達関数の算出

G_z1w = 1/s*S*Ps;
G_z2w = 0.1*S*Ps*Ks;

G_zw = [G_z1w G_z2w];
%% 
% 特異値プロット

sigma(G_zw)
%% 
% お手製特異値プロット

% 周波数範囲の定義（例：0 から 100 rad/s）
omega = logspace(-1, 2, 500);

% 各周波数での特異値を格納するための配列
singular_values = zeros(1, length(omega));

% 各周波数における特異値の計算
for i = 1:length(omega)
    % 伝達関数行列の周波数応答を計算
    [mag,~,~] = bode(G_zw, omega(i));
    
    % 特異値の計算
    singular_values(i) = max(svd(squeeze(mag)));
end

Hinf_norm = max(singular_values)

semilogx(omega,20*log10(singular_values))
%% 
%