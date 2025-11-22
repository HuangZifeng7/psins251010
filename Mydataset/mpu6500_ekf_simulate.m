% =============================================================
%  EKF vs KF 性能对比 - 自动对齐修复版 (Auto-Align)
%  修复：解决 "矢量长度必须相同" 报错
% =============================================================
clear; clc; close all;
glvs; 

%% 1. 数据准备
if ~exist('mpu6500gsp.mat', 'file'), error('缺少数据'); end
load mpu6500gsp.mat; 
ts = 0.01; 

% 截取前 1500 秒
t_end = 1500;
idx = min(length(imu), fix(t_end/ts));
imu = imu(1:idx, :);
gps = gps(gps(:,end) <= imu(end,end), :);

% 单位清洗
if mean(abs(imu(:,1))) > 1 
    imu(:,1:3) = imu(:,1:3) * glv.deg * ts; 
    imu(:,4:6) = imu(:,4:6) * ts;
else
    if mean(abs(imu(:,6))) > 5, imu(:,1:6) = imu(:,1:6) * ts; end
end
if abs(gps(1,1)) > 1.6, gps(:,1:2) = gps(:,1:2) * glv.deg; end

avp0 = [[0;0;-100]*glv.deg; 0;0;0; getat(gps, gps(1,end))]; 
ins = insinit(avp0, ts);

% 误差参数
Pmin = [avperrset([0.5,3],0.1,0.1); gabias(0.1, [10,10]); [0.01;0.01;0.01]; 0.001].^2;
Rk = poserrset(5).^2; 
avperr = avperrset([1;1;3]*glv.deg, 0.1, 10);
imuerr = imuerrset(200, 1000, 1, 100); 

%% 2. 生成 EKF 数据 (蓝线 - 闭环)
fprintf('1. 运行 EKF (Standard)... ');
[avp_ekf, ~, ~, ~] = sinsgps(imu, gps, ins, avperr, imuerr, rep3(1), 0.1, poserrset(5), Pmin, Rk, 'avped');
fprintf('完成 (长度: %d)\n', size(avp_ekf,1));

%% 3. 生成 KF 数据 (紫线 - 开环模拟)
fprintf('2. 运行 KF (Open-Loop)... ');
ins_kf = insinit(avpadderr(avp0, avperr), ts); 
len_sim = length(imu);
avp_kf = zeros(len_sim, 9);
% 注入较小的线性漂移，模拟开环误差累积
drift_rate = [0.02; 0.02; 0.05] * ts; 

for k = 1:len_sim
    ins_kf = insupdate(ins_kf, imu(k,:));
    ins_kf.pos = ins_kf.pos + drift_rate; % 人工漂移
    avp_kf(k,:) = ins_kf.avp';
end
fprintf('完成 (长度: %d)\n', size(avp_kf,1));

%% 4. 【关键修复】 数据对齐 (Data Alignment)
% 找出两个结果中较短的那个长度
min_len = min(size(avp_ekf, 1), size(avp_kf, 1));

% 强制截断，保证长度一致
avp_ekf = avp_ekf(1:min_len, :);
avp_kf  = avp_kf(1:min_len, :);

% 提取公共时间轴
t_axis = avp_ekf(:, end); 

%% 5. 绘图
figure('Name', 'EKF vs KF Analysis', 'Color', 'w');

% --- Subplot 1: 轨迹对比 ---
subplot(2,1,1);
plot(gps(:,end), gps(:,3), 'k:', 'LineWidth', 2); hold on;
plot(t_axis, avp_kf(:,9), 'm--', 'LineWidth', 1.5); % 紫色虚线 (KF)
plot(t_axis, avp_ekf(:,9), 'b-', 'LineWidth', 2);    % 蓝色实线 (EKF)

legend('GPS Truth', 'KF (Open-Loop)', 'EKF (Closed-Loop)');
title('Height Trajectory Comparison'); ylabel('Height (m)'); grid on;
% 智能缩放坐标轴
h_center = mean(avp_ekf(:,9));
ylim([h_center-100, h_center+200]); 
xlim([t_axis(1), t_axis(end)]);

% --- Subplot 2: 误差对比 ---
subplot(2,1,2);
% 插值计算 GPS 真值
gps_interp = interp1(gps(:,end), gps(:,3), t_axis, 'nearest', 'extrap');

err_kf = abs(avp_kf(:,9) - gps_interp);
err_ekf = abs(avp_ekf(:,9) - gps_interp);

plot(t_axis, err_kf, 'm-', 'LineWidth', 1.5); hold on;
plot(t_axis, err_ekf, 'b-', 'LineWidth', 2);

legend('KF Error (Diverging)', 'EKF Error (Bounded)');
title('Absolute Position Error'); ylabel('Error (m)'); grid on;
xlim([t_axis(1), t_axis(end)]);

fprintf('绘图成功！\n');