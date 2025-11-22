% =============================================================
%  第一层 KF 精度验证 - 量化分析版
%  功能：对比 Raw vs Filtered 的 STD(噪声) 和 漂移误差
% =============================================================
clear; clc; close all;
glvs; 

%% 1. 数据加载
if ~exist('mpu6500gsp.mat', 'file'), error('缺少数据'); end
load mpu6500gsp.mat; 
ts = 0.01; 

% 截取前 200秒 (验证滤波效果不需要太长，短时漂移更能看出噪声影响)
t_end = 200; 
idx = min(length(imu), fix(t_end/ts));
imu_raw = imu(1:idx, :);
gps = gps(gps(:,end) <= imu_raw(end,end), :);

% 单位清洗
if mean(abs(imu_raw(:,1))) > 1 
    imu_raw(:,1:3) = imu_raw(:,1:3) * glv.deg * ts; 
    imu_raw(:,4:6) = imu_raw(:,4:6) * ts;
end
if abs(gps(1,1)) > 1.6, gps(:,1:2) = gps(:,1:2) * glv.deg; end

%% 2. 第一层 KF 滤波 (传感器级)
fprintf('正在执行 KF 滤波...\n');
imu_clean = zeros(size(imu_raw));
len = length(imu_raw);

% 参数设置 (根据上一轮的经验)
Q_gyro = 1e-6;  R_gyro = 0.1;
Q_acc  = 1e-5;  R_acc  = 0.5;

for axis = 1:6
    if axis <= 3, Q=Q_gyro; R=R_gyro; else, Q=Q_acc; R=R_acc; end
    
    z_in = imu_raw(:, axis);
    x_out = zeros(len, 1);
    x = z_in(1); P = R;
    
    for k = 1:len
        % 标量 KF
        x_pred = x; P_pred = P + Q;
        z = z_in(k);
        K = P_pred / (P_pred + R);
        x = x_pred + K * (z - x_pred);
        P = (1 - K) * P_pred;
        x_out(k) = x;
    end
    imu_clean(:, axis) = x_out;
end

%% 3. 【核心】量化指标计算 1：噪声抑制能力 (STD)
% 我们选取一段相对平稳的数据计算 STD (假设前100个点是静止或平稳的)
N_static = 100;

% 陀螺仪 X 轴
std_gyro_raw = std(imu_raw(1:N_static, 1));
std_gyro_cln = std(imu_clean(1:N_static, 1));
improve_gyro = (std_gyro_raw - std_gyro_cln) / std_gyro_raw * 100;

% 加速度计 Z 轴
std_acc_raw = std(imu_raw(1:N_static, 6));
std_acc_cln = std(imu_clean(1:N_static, 6));
improve_acc = (std_acc_raw - std_acc_cln) / std_acc_raw * 100;

fprintf('\n=== 1. 信号层精度对比 (Signal Precision) ===\n');
fprintf('陀螺仪 STD: 原始 = %.5f -> 滤波后 = %.5f (降低了 %.2f%%)\n', std_gyro_raw, std_gyro_cln, improve_gyro);
fprintf('加速度 STD: 原始 = %.5f -> 滤波后 = %.5f (降低了 %.2f%%)\n', std_acc_raw, std_acc_cln, improve_acc);

%% 4. 纯惯导解算 (为了验证轨迹精度)
fprintf('\n正在进行纯惯导解算对比...\n');

% 扣除零偏 (为了公平对比噪声影响)
bias_g = mean(imu_raw(1:N_static, 1:3));
bias_a = mean(imu_raw(1:N_static, 4:6)); bias_a(3) = bias_a(3) - glv.g0*ts;

avp0 = [[0;0;-100]*glv.deg; 0;0;0; getat(gps, gps(1,end))]; 
ins_raw = insinit(avp0, ts);
ins_cln = insinit(avp0, ts);

avp_res_raw = zeros(len, 9);
avp_res_cln = zeros(len, 9);

for k = 1:len
    % 原始数据 (扣零偏)
    dat_r = imu_raw(k,:);
    dat_r(1:3) = dat_r(1:3) - bias_g; dat_r(4:6) = dat_r(4:6) - bias_a;
    ins_raw = insupdate(ins_raw, dat_r);
    avp_res_raw(k,:) = ins_raw.avp';
    
    % 滤波数据 (扣零偏)
    dat_c = imu_clean(k,:);
    dat_c(1:3) = dat_c(1:3) - bias_g; dat_c(4:6) = dat_c(4:6) - bias_a;
    ins_cln = insupdate(ins_cln, dat_c);
    avp_res_cln(k,:) = ins_cln.avp';
end

%% 5. 【核心】量化指标计算 2：位置精度 (RMSE / Max Error)
t = (1:len)*ts;
gps_interp = interp1(gps(:,end), gps(:,3), t, 'nearest', 'extrap');

% 计算绝对误差
err_raw = abs(avp_res_raw(:,9) - gps_interp');
err_cln = abs(avp_res_cln(:,9) - gps_interp');

% 计算 RMSE (均方根误差)
rmse_raw = sqrt(mean(err_raw.^2));
rmse_cln = sqrt(mean(err_cln.^2));
improve_rmse = (rmse_raw - rmse_cln) / rmse_raw * 100;

% 计算最大漂移 (Max Drift)
max_raw = max(err_raw);
max_cln = max(err_cln);

fprintf('\n=== 2. 轨迹层精度对比 (Trajectory Precision) ===\n');
fprintf('高度 RMSE: 原始 = %.2fm -> 滤波后 = %.2fm (精度提升 %.2f%%)\n', rmse_raw, rmse_cln, improve_rmse);
fprintf('最大漂移 : 原始 = %.2fm -> 滤波后 = %.2fm\n', max_raw, max_cln);

%% 6. 绘图
figure('Name', 'Verification of 1st Layer KF', 'Color', 'w');

% 图1: 信号对比
subplot(2,2,1);
plot(t(1:500), imu_raw(1:500,1)/ts/glv.deg, 'Color', [0.7 0.7 0.7]); hold on;
plot(t(1:500), imu_clean(1:500,1)/ts/glv.deg, 'b', 'LineWidth',1.5);
title('Signal Noise Reduction (Gyro X)'); ylabel('deg/s'); grid on;

% 图2: 轨迹对比
subplot(2,2,[2,4]); % 右侧大图
plot(gps(:,end), gps(:,3), 'k:', 'LineWidth', 2); hold on;
plot(t, avp_res_raw(:,9), 'm--', 'LineWidth', 1.5);
plot(t, avp_res_cln(:,9), 'b-', 'LineWidth', 2);
legend('Ground Truth', 'INS (Raw)', 'INS (Filtered)');
title(['Trajectory Drift (RMSE Improved by ' num2str(improve_rmse, '%.1f') '%)']); 
ylabel('Height (m)'); grid on;

% 图3: 误差对比
subplot(2,2,3);
plot(t, err_raw, 'm', 'LineWidth', 1); hold on;
plot(t, err_cln, 'b', 'LineWidth', 2);
title('Position Error'); ylabel('Error (m)'); legend('Raw Error', 'Filtered Error'); grid on;