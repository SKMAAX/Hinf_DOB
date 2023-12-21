%% 拡張H∞制御のLMI解法で外乱オブザーバ設計
% バネマスダンパ(Mass-Spring-Damper)系への適用

clear;
s = tf('s');
%% 
% 相補感度関数Taに対する重みWaの設定

% Wa = tf([0.1 0.1],[0.01 1]);
theta = 10;
Aa = -100;
ba = -990;
ca = 1;
[na,~] = size(Aa);

Wa = ss(Aa, ba, ca, theta);
bodemag(Wa); grid on;
%% 
% 制御対象の記述

a0 = 1; % ばね要素
a1 = 0.5; % ダンパ要素
A0 = [0 1; -a0 -a1];
b = [0;1]; 
c = [1 0]; % 位置観測
[n0,~] = size(A0);

%% 
% 一般化制御対象の記述

A = [A0 b zeros(n0,1);zeros(1,n0) 0 0; zeros(1,n0) 0 Aa];
[n,~] = size(A);
B1 = [zeros(n0,2);1 0;0 ba];
B2 = [b;0;0];
C1 = [zeros(1,n0) 1 0];
C2 = [c 0 ca];
D12 = [1];
D21 = [0 theta];
D12_dag = 1;
D21_dag = [0;1/theta];

P = ss(A,B2,C2,0);
Ps = tf(P(1));
bode(Ps,Ps+Wa); %制御対象と、制御対象＋不確かさワーストケースのボード線図
%% 
% LMIで使用する各行列の定義

U0 = zeros(n,n);
U1 = eye(n);
U1_dag = pinv(U1);
A_U = U1*A*U1_dag;
C1_U = C1*U1_dag;
C2_U = C2*U1_dag;
B1_U = U1*B1;
U_L = zeros(2,1);
%% 
% オブザーバゲインの算出

% SDPソルバとLMIパーサのパスの設定
addpath(genpath('～～sedumi-masterへのパス～～'))
addpath(genpath('～～YALMIP-masterへのパス～～'))

Y_U=sdpvar(n,n,'sy');
N=sdpvar(n,1); % N=Y_U*U_1*L;
gamma = sdpvar(1);

% -------------------
LMI=[]; % LMI initialized

ep = 1e-6;
LMI=[LMI, Y_U>=ep*eye(n)];

LMI2 = [A_U'*Y_U+Y_U*A_U+C2_U'*N'+N*C2_U+C1_U'*C1_U   Y_U*B1_U+N*D21;
    B1_U'*Y_U+D21'*N'  -eye(2)*gamma];
LMI=[LMI, LMI2<=-ep*eye(n+2)];

% -------------------
% ops=sdpsettings; ops.shift=1e-6;
sol=solvesdp(LMI,gamma);

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

K = ss(A+B2*F+L*C2,[L B2],[-F;C2],0);
tfK = tf(K);
Ks = tfK(1,1);
%% 
% 感度関数と相補感度関数

S = tf(feedback(1,Ks*Ps,+1)); %inv(1-P*K);
T = tf(feedback(Ks*Ps,1,+1)); %inv(1-P*K)*P*K;
Ta = tf(feedback(Ks,Ps,+1)); %inv(1-P*K)*K;
%% 
% 準相補感度関数Taと、重みWaの大小関係確認

bodemag(Ta,1/Wa); grid on;
title('準相補感度関数Taと、重みWaの大小関係確認')
legend('Ta','1/Wa')
%% 
% ロバスト安定性の限界について確認

margin(1.1*Ta*Wa) % …ゲイン余裕&位相余裕はないが、ゲイン1.1倍で即破綻するわけでは無さそう。

delays = tf(1,1,'inputDelay',0.001); % Ta*Waが位相線図的にまだ余裕ありそうなので、適当なディレイを挿入
bode(Ta*Wa,1.1*Ta*Wa*delays,{0.1,10000})
title('想定内ワーストケースでの不確かさ込み閉ループと、想定を超えて与えた例')
legend('Ta*Wa','1.1*Ta*Wa*delay','Location','southwest')

% Simulinkモデルの実行

simumodel = 'DOB_MSD_LMI';
open_system(simumodel)
set_param(simumodel+"/不確かさ切替設定",'Value','0') %不確かさなし設定
set_param(simumodel+"/制御入力有無",'sw','0') %制御入力なし設定
out_noctrl = sim(simumodel);
%% 
% 制御入力なしでの結果表示

lw = 1;
FontSize = 10;
subplot(2,1,1),plot(out_noctrl.d.Time,out_noctrl.d.Data,out_noctrl.est_d.Time,out_noctrl.est_d.Data,...
    'LineWidth',lw); grid on;
leg1 = legend('d', 'estimated d');
leg1.FontSize = FontSize;
set(leg1,'Location','SouthEast');
title('w/o Feedback, w/o Uncertainty','FontSize',FontSize)
ylabel('[N]','FontSize',FontSize)
xlabel('Time[s]','FontSize',FontSize)
%% 
% 条件の切り替え、Simulinkモデルの実行

set_param(simumodel+"/制御入力有無",'sw','1') %制御入力あり設定
set_param(simumodel+"/不確かさ切替設定",'Value','0') %不確かさなし設定
out_ctrl = sim(simumodel);
%% 
% 制御入力ありでの結果表示 → 入力なしと同一の結果

subplot(2,1,2),plot(out_noctrl.d.Time,out_noctrl.d.Data,out_noctrl.est_d.Time,out_noctrl.est_d.Data,...
    'LineWidth',lw); grid on;
leg2 = legend('d', 'estimated d');
leg2.FontSize = FontSize;
set(leg2,'Location','SouthEast');
title('w/ Feedback, w/o Uncertainty','FontSize',FontSize)
ylabel('[N]','FontSize',FontSize)
xlabel('Time[s]','FontSize',FontSize)
%% 
% 条件の切り替え、Simulinkモデルの実行

set_param(simumodel+"/制御入力有無",'sw','0') %制御入力なし設定
set_param(simumodel+"/不確かさ切替設定",'Value','1') %不確かさあり設定
out_uncert = sim(simumodel);
%% 
% 不確かさあり、制御入力なしでの結果表示 

figure;
subplot(2,1,1),plot(out_uncert.d.Time,out_uncert.d.Data,out_uncert.est_d.Time,out_uncert.est_d.Data,...
    'LineWidth',lw); grid on;
leg3 = legend('d', 'estimated d');
leg3.FontSize = FontSize;
set(leg3,'Location','SouthEast');
title('w/o Feedback, w/ Uncertainty','FontSize',FontSize)
ylabel('[N]','FontSize',FontSize)
xlabel('Time[s]','FontSize',FontSize)
ylim([0,1.4])
%% 
% 条件の切り替え、Simulinkモデルの実行

set_param(simumodel+"/制御入力有無",'sw','1') %制御入力あり設定
set_param(simumodel+"/不確かさ切替設定",'Value','1') %不確かさあり設定
out_ctrl_uncert = sim(simumodel);
%% 
% 不確かさあり、制御入力ありでの結果表示 

subplot(2,1,2),plot(out_ctrl_uncert.d.Time,out_ctrl_uncert.d.Data,out_ctrl_uncert.est_d.Time,out_ctrl_uncert.est_d.Data,...
    'LineWidth',lw); grid on;
leg4 = legend('d', 'estimated d');
leg4.FontSize = FontSize;
set(leg4,'Location','SouthEast');
title('w/ Feedback, w/ Uncertainty','FontSize',FontSize)
ylabel('[N]','FontSize',FontSize)
xlabel('Time[s]','FontSize',FontSize)
ylim([0,1.4])
%% 
% 想定の範囲を超えた不確かさに対する挙動の確認

set_param(simumodel+"/制御入力有無",'sw','0') %制御入力なし設定
set_param(simumodel+"/不確かさ切替設定",'Value','2') %想定以上の不確かさあり設定
set_param(simumodel+"/不確かさゲイン",'Gain','1.1') %不確かさ許容オーバー量
set_param(simumodel+"/遅延",'DelayTime','0.001') %不確かさ許容オーバー量
out_overuncert = sim(simumodel);
%% 
% 想定を超える不確かさ、制御入力なしでの結果表示 

figure;
subplot(2,1,1),plot(out_overuncert.d.Time,out_overuncert.d.Data,out_overuncert.est_d.Time,out_overuncert.est_d.Data,...
    'LineWidth',lw); grid on;
leg5 = legend('d', 'estimated d');
leg5.FontSize = FontSize;
set(leg5,'Location','SouthEast');
title('w/o Feedback, w/ Large Uncertainty','FontSize',FontSize)
ylabel('[N]','FontSize',FontSize)
xlabel('Time[s]','FontSize',FontSize)
ylim([0,1.4])

%% 
% 条件の切り替え、Simulinkモデルの実行

set_param(simumodel+"/制御入力有無",'sw','1') %制御入力あり設定
out_ctrl_overuncert = sim(simumodel);
%% 
% 想定を超える不確かさ、制御入力ありでの結果表示 

subplot(2,1,2),plot(out_ctrl_overuncert.d.Time,out_ctrl_overuncert.d.Data,out_ctrl_overuncert.est_d.Time,out_ctrl_overuncert.est_d.Data,...
    'LineWidth',lw); grid on;
leg6 = legend('d', 'estimated d');
leg6.FontSize = FontSize;
set(leg6,'Location','SouthEast');
title('w/ Feedback, w/ Large Uncertainty','FontSize',FontSize)
ylabel('[N]','FontSize',FontSize)
xlabel('Time[s]','FontSize',FontSize)
ylim([0,1.4])