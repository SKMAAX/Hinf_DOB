%% 極配置法による外乱オブザーバ設計
% バネマスダンパ(Mass-Spring-Damper)系への適用

clear;
s = tf('s');
%% 
% 相補感度関数Taに対する重みWaの設定

%%不確かさパターン1
% Wa = tf([0.1 0.1],[0.01 1]); 
theta = 10;
Aa = -100;
Ba = -990;
Ca = 1;
%%不確かさパターン2
% Aa = [-16.05  -15.94;16.08   15.91];
% Ba = [1.52;-1.52];
% Ca = [1.747  1.73]*5;
% theta = 0.00578*5;
[na,~] = size(Aa);

Wa = ss(Aa, Ba, Ca, theta);
bodemag(Wa);
%% 
% 制御対象の記述

a0 = 1; % ばね要素
a1 = 0.5; % ダンパ要素
A0 = [0 1; -a0 -a1];
b = [0;1]; 
c = [1 0]; % 位置観測
[n0,~] = size(A0);

P = ss(A0,b,c,0);
%% 
% 一定値外乱を含めた拡大系の記述

A = [A0 b ;zeros(1,n0) 0];
B = [b; 0];
C = [c 0];


%% 
% オブザーバゲインの算出

p = -[20; 22; 24];

L = -place(A',C',p)';
%% 
% 制御器K(s)の算出

% K = ss(A+B2*F+L*C2,[L B2],[-F;C2],[0 -1;1 0]);
% tfK = tf(K);
% Ks = tfK(1,1);
%% 
% 感度関数と相補感度関数

% S = tf(feedback(1,Ks*Ps,+1)); %inv(1-P*K);
% Ta = tf(feedback(Ks,Ps,+1)); %inv(1-P*K)*K;
%% 
% 相補感度関数Tと、重み関数Waの大小関係確認

% bodemag(Ta,1/Wa);
% Simulinkモデルの実行

simumodel = 'DOB_MSD_poleplace';
open_system(simumodel)
set_param(simumodel+"/不確かさ(ワーストケース)の有無",'sw','0') %不確かさなし設定
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
set_param(simumodel+"/不確かさ(ワーストケース)の有無",'sw','0') %不確かさなし設定
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
set_param(simumodel+"/不確かさ(ワーストケース)の有無",'sw','1') %不確かさあり設定
out_uncert = sim(simumodel);
%% 
% 不確かさあり、制御入力なしでの結果表示 → 発散する

figure;
subplot(2,1,1),plot(out_uncert.d.Time,out_uncert.d.Data,out_uncert.est_d.Time,out_uncert.est_d.Data,...
    'LineWidth',lw); grid on;
leg3 = legend('d', 'estimated d');
leg3.FontSize = FontSize;
set(leg3,'Location','NorthEast');
title('w/o Feedback, w/ Uncertainty','FontSize',FontSize)
ylabel('[N]','FontSize',FontSize)
xlabel('Time[s]','FontSize',FontSize)
ylim([0,5])
%% 
% 条件の切り替え、Simulinkモデルの実行

set_param(simumodel+"/制御入力有無",'sw','1') %制御入力あり設定
set_param(simumodel+"/不確かさ(ワーストケース)の有無",'sw','1') %不確かさあり設定
out_ctrl_uncert = sim(simumodel);
%% 
% 不確かさあり、制御入力ありでの結果表示 → 発散する

subplot(2,1,2),plot(out_ctrl_uncert.d.Time,out_ctrl_uncert.d.Data,out_ctrl_uncert.est_d.Time,out_ctrl_uncert.est_d.Data,...
    'LineWidth',lw); grid on;
leg4 = legend('d', 'estimated d');
leg4.FontSize = FontSize;
set(leg4,'Location','NorthEast');
title('w/ Feedback, w/ Uncertainty','FontSize',FontSize)
ylabel('[N]','FontSize',FontSize)
xlabel('Time[s]','FontSize',FontSize)
ylim([0,5])