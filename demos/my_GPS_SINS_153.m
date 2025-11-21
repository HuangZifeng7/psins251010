%% WHU-i2Nav GINS 数据集 SINS/GNSS 组合导航解算脚本
% 适配数据: whu_gins_data.mat (由之前的脚本转换得到)
% 核心逻辑: 15维误差状态卡尔曼滤波 (松耦合)
% Time: 2025
clc; clear; close all;
glvs; 
psinstypedef(153); % 使用15维状态模型

%% 1. 加载数据
disp('正在加载数据...');
if ~exist('whu_gins_data_20251121_171945.mat', 'file')
    error('找不到 whu_gins_data.mat，请先运行数据转换脚本！');
end
load('whu_gins_data.mat'); 
% 变量: imu(FRD delta-angle/delta-velocity, time), avp_ref(真值), gps_data(观测), ts

% --- 诊断：检查 imu_type 与时间范围 ---
if ~exist('imu','var')
    error('变量 imu 未在 whu_gins_data.mat 中找到，请先运行数据转换脚本');
end
if ~exist('imu_type','var')
    imu_type = 'unknown';
    warning('未在数据文件中找到 imu_type 字段，默认设为 unknown');
end
disp(['载入数据完成：imu_type = ', imu_type]);
gyr_mean = mean(abs(imu(:,1:3)),1);
acc_mean = mean(abs(imu(:,4:6)),1);
disp(['IMU 陀螺平均绝对值 (per-axis): ', num2str(gyr_mean)]);
disp(['IMU 加速度平均绝对值 (per-axis): ', num2str(acc_mean)]);
disp(['IMU 时间范围: ', num2str(imu(1,end)), ' -> ', num2str(imu(end,end))]);
disp(['GPS 时间范围: ', num2str(gps_data(1,end)), ' -> ', num2str(gps_data(end,end))]);

%% 2. 初始化设置
disp('正在初始化...');

% [重要] 截取时间段：为了演示，我们从第1秒跑到最后
% 你可以通过修改 start_t 和 end_t 来只跑一部分数据
start_t = imu(1, end); 
end_t   = imu(end, end); 

% 找到对应的索引范围
idx_s = find(imu(:,end) >= start_t, 1);
idx_e = find(imu(:,end) <= end_t, 1, 'last');

% 截取数据
imu_run = imu(idx_s:idx_e, :);

% === 初始化 INS 状态 ===
% 直接使用真值(avp_ref)的第一帧作为初始状态
% 实际工程中需要进行粗对准(align)，但这里为了验证算法，先用真值启动
% avp_ref 格式: [Att(3), Vel(3), Pos(3), Time]
p0_ref = avp_ref(idx_s, 1:9)'; 
ins = insinit(p0_ref, ts);

% === [非常重要] 杆臂值设置 (Lever Arm) ===
% 根据 README，ADIS16465 的杆臂 (FRD frame)
% 如果不加这个，组合导航的位置误差会有几十厘米的系统偏差
ins.ant = [-0.073; 0.302; 0.087]; 

% === 卡尔曼滤波初始化 ===
% 根据 README 中的 ADIS16465 参数设置噪声
% Gyro Bias: 25 deg/h, Acc Bias: 200 mGal
% ARW: 0.1 deg/sqr(h), VRW: 0.1 m/s/sqr(h)
imuerr = imuerrset(25, 200, 0.1, 0.1); 

% 初始误差协方差 P0 (假设初始给定的位置速度比较准，姿态有点误差)
davp0 = avperrset([1;1;5], 0.1, [1;1;1]); 

% GNSS 量测噪声 R (RTK 精度很高，设小一点)
rk = poserrset([0.05; 0.05; 0.1]); 

kf = kfinit(ins, davp0, imuerr, rk);

% 过程噪声 Q 阵设置 (根据 IMU 噪声水平)
kf.Pmin = [avperrset(0.01,1e-4,0.1); gabias(1e-3, [1,10])].^2;
kf.Qk = diag([imuerr.web; imuerr.wdb; zeros(9,1)])^2 * ts; 

%% 3. 循环解算 (Filter Loop)
len = length(imu_run);
[nn, ts, nts] = nnts(2, ts); % 2子样算法
k1 = 1; % 结果数组索引
gps_ptr = 1; % GPS 数据指针

% 预分配内存
res_avp = zeros(fix(len/nn), 10); % 存储解算结果 [Att, Vel, Pos, T]
res_xk  = zeros(fix(len/nn), 16); % 存储误差状态 [Phi, dV, dP, T]
% 预分配 GPS-update 日志（每次 GPS 更新的量测/残差）
gpsN = size(gps_data,1);
allZ_rad = nan(gpsN,3);    % Z in rad/mix (dLat rad, dLon rad, dH m)
allZ_m   = nan(gpsN,3);    % converted to meters [dN,dE,dH]
allRk    = nan(gpsN,max(3,kf.m)); % residuals recorded (kf.m should be <=3 for pos)
allKnorm = nan(gpsN,1);    % norm of Kalman gain
gps_log_idx = 0;

timebar(nn, len, 'SINS/GNSS Processing...');

for k = 1:nn:len-nn+1
    % --- 1. 惯导机械编排 (INS Mechanization) ---
    k_end = k+nn-1;
    wvm = imu_run(k:k_end, 1:6); % 提取增量
    t_now = imu_run(k_end, end); % 当前时间
    
    ins = insupdate(ins, wvm);
    
    % --- 2. KF 时间更新 (Prediction) ---
    kf.Phikk_1 = kffk(ins);
    kf = kfupdate(kf);
    
    % --- 3. KF 量测更新 (Correction) ---
    % 检查是否有 GPS 数据到来
    % 逻辑：如果 GPS 指针没超限，且当前 IMU 时间 超过了 GPS 数据的时间戳
    % 支持在同一 IMU 步内到达多条 GPS 的情况
    while gps_ptr <= size(gps_data, 1) && t_now >= gps_data(gps_ptr, end)
        % 获取当前时刻的 GPS 观测值
        pos_gps = gps_data(gps_ptr, 4:6)'; % Lat(rad), Lon(rad), Hgt(m)
        vn_gps  = gps_data(gps_ptr, 1:3)'; % Vn, Ve, Vd

        % 计算 INS-GPS 差并打印（以米为单位，便于直观判断）
        Z = ins.pos - pos_gps; % [dLat(rad); dLon(rad); dH(m)]
        latmean = (ins.pos(1) + pos_gps(1))/2;
        dN = Z(1) * glv.Re;                      % 北向差, m
        dE = Z(2) * glv.Re * cos(latmean);      % 东向差, m
        dH = Z(3);                               % 高度差, m
        %disp(['正在执行第 ', num2str(gps_ptr), ' 次 GPS 更新: ΔN=', sprintf('%.3f',dN), ' m, ΔE=', sprintf('%.3f',dE), ' m, ΔH=', sprintf('%.3f',dH), ' m']);

        % 记录并执行量测更新 (这里只用了位置更新)
        gps_log_idx = gps_log_idx + 1;
        allZ_rad(gps_log_idx,:) = Z';
        % convert to meters for logging
        allZ_m(gps_log_idx,1) = dN; allZ_m(gps_log_idx,2) = dE; allZ_m(gps_log_idx,3) = dH;

        kf = kfupdate(kf, Z, 'M');
        % 记录卡尔曼残差和 K 矩阵指标（若存在）
        if isfield(kf,'rk') && ~isempty(kf.rk)
            r = kf.rk(:)'; allRk(gps_log_idx,1:numel(r)) = r; end
        if isfield(kf,'Kk') && ~isempty(kf.Kk)
            allKnorm(gps_log_idx) = norm(kf.Kk,'fro'); end

        [kf, ins] = kffeedback(kf, ins, 1, 'avp');

        gps_ptr = gps_ptr + 1;
    end
    
    % --- 4. 保存结果 ---
    res_avp(k1, :) = [ins.avp', t_now];
    res_xk(k1, :)  = [kf.xk(1:15)', t_now]; 
    k1 = k1 + 1;
    
    timebar;
end
res_avp(k1:end, :) = []; % 删除多余空间
res_xk(k1:end, :)  = [];

disp('解算完成！');

%% 4. 绘图与误差分析
% 对齐真值 (Resample reference to result time)
% 因为解算结果的频率可能跟真值不完全一样，用插值对齐
% --- 插值真值到解算结果时间点（更稳健） ---
% 首先确保 avp_ref 时间列单调（之前已检查）并去除任何包含 NaN 的源行
valid_idx = all(~isnan(avp_ref(:,1:9)),2);
if ~all(valid_idx)
    warning('avp_ref 包含 NaN 行，已丢弃 %d 行用于插值', sum(~valid_idx));
end
T_ref = avp_ref(valid_idx,end);
V_ref = avp_ref(valid_idx,1:9);

% 使用线性插值并允许外推 (extrap) 以避免未覆盖时间段导致的 NaN
avp_ref_interp = interp1(T_ref, V_ref, res_avp(:,end), 'linear', 'extrap');

% 如果插值后仍包含 NaN（非常少见），用最近邻填充
if any(isnan(avp_ref_interp(:)))
    warning('插值后仍含 NaN，应用最近邻填充');
    % 对每列进行 fillmissing with 'nearest'
    for c=1:size(avp_ref_interp,2)
        avp_ref_interp(:,c) = fillmissing(avp_ref_interp(:,c),'nearest');
    end
end

% 调用 PSINS 绘图工具
% 1. 轨迹对比
% --- 诊断：打印样例数据并绘制 2D 轨迹（使用 pos 列） ---
disp('--- 数据样例 (首5行) ---');
disp('avp_ref(1:5,:)'); disp(avp_ref(1:5,:));
disp('gps_data(1:5,:)'); disp(gps_data(1:5,:));
disp('p0_ref (initial avp used to init ins):'); disp(p0_ref');

figure;
% 使用 gps_data 的 Lon/Lat (gps_data: [Vn Ve Vd Lat Lon Hgt Time])
plot(gps_data(:,5)/glv.deg, gps_data(:,4)/glv.deg, 'k.'); hold on; % GNSS observed (Lon,Lat)
% 使用 res_avp 的 pos 列 (res_avp: [att(1:3), vel(1:3), pos(1:3), time])
plot(res_avp(:,8)/glv.deg, res_avp(:,7)/glv.deg, 'r');
legend('GNSS Observed', 'SINS/GNSS Integrated');
xlabel('Longitude (deg)'); ylabel('Latitude (deg)'); grid on;
title('2D Trajectory');

% 2. 误差曲线
avp_err = res_avp(:, 1:9) - avp_ref_interp;
% 姿态误差转为角分 (arcmin)
avp_err(:,1:3) = avp_err(:,1:3) / glv.min; % arcmin
% 速度误差保持 m/s，位置差将转换为米用于显示

figure;
subplot(3,1,1); plot(res_avp(:,end), avp_err(:,1:3)); grid on;
title('Attitude Error (arcmin)'); legend('\phi', '\theta', '\psi');
subplot(3,1,2); plot(res_avp(:,end), avp_err(:,4:6)); grid on;
title('Velocity Error (m/s)'); legend('Vn', 'Ve', 'Vd');
% 将经纬弧度差转换为米
meanlat = mean(avp_ref_interp(:,7));
pos_err_rad = avp_err(:,7:9); % [dLat(rad), dLon(rad), dH(m)]
pos_err_m = zeros(size(pos_err_rad));
pos_err_m(:,1) = pos_err_rad(:,1) * glv.Re;                % dN (m)
pos_err_m(:,2) = pos_err_rad(:,2) * glv.Re * cos(meanlat); % dE (m)
pos_err_m(:,3) = pos_err_rad(:,3);                         % dD (m)
subplot(3,1,3); plot(res_avp(:,end), pos_err_m); grid on;
title('Position Error (m)'); legend('dN', 'dE', 'dD');
% 输出统计
% 计算 RMS 时忽略 NaN（更健壮）
rmsN = sqrt(nanmean(pos_err_m(:,1).^2));
rmsE = sqrt(nanmean(pos_err_m(:,2).^2));
rmsH = sqrt(nanmean(pos_err_m(:,3).^2));
fprintf('\nPosition RMS (m) (nan-ignoring): N=%.3f, E=%.3f, H=%.3f\n', rmsN, rmsE, rmsH);

% --- 输出并绘制 GPS 更新 residuals 与量测统计 ---
if gps_log_idx>0
    allZ_m = allZ_m(1:gps_log_idx,:);
    allZ_rad = allZ_rad(1:gps_log_idx,:);
    allRk = allRk(1:gps_log_idx,1:3);
    allKnorm = allKnorm(1:gps_log_idx);

    fprintf('GPS updates logged: %d\n', gps_log_idx);
    fprintf('Mean |Z| (m): N=%.3f E=%.3f H=%.3f\n', mean(abs(allZ_m(:,1))), mean(abs(allZ_m(:,2))), mean(abs(allZ_m(:,3))));
    fprintf('RMS residuals (m for pos-components rad->m): N=%.3f E=%.3f H=%.3f\n', sqrt(nanmean((allRk(:,1)*glv.Re).^2)), sqrt(nanmean((allRk(:,2)*glv.Re*cos(mean(avp_ref(:,7)))).^2)), sqrt(nanmean((allRk(:,3)).^2)));

    figure; subplot(3,1,1); histogram(allZ_m(:,1),50); title('GPS Z dN (m) histogram'); subplot(3,1,2); histogram(allZ_m(:,2),50); title('GPS Z dE (m) histogram'); subplot(3,1,3); histogram(allZ_m(:,3),50); title('GPS Z dH (m) histogram');
    figure; subplot(2,1,1); histogram(allRk(:,1)*glv.Re,50); title('Residual rk dN (m) histogram'); subplot(2,1,2); histogram(allRk(:,2)*glv.Re*cos(mean(avp_ref(:,7))),50); title('Residual rk dE (m) histogram');
    figure; plot(allKnorm); title('KF gain norm per GPS update'); grid on;

    % 保存 debug 信息
    out_debug = ['out_debug_', datestr(now,'yyyymmdd_HHMMSS'), '.mat'];
    save(out_debug, 'allZ_rad','allZ_m','allRk','allKnorm','gps_log_idx');
    disp(['已保存 GPS 更新诊断数据到: ', out_debug]);
end

% 3. 陀螺仪零偏估计曲线

figure;
plot(res_xk(:,end), res_xk(:,10:12)/glv.dph); grid on;
title('Estimated Gyro Bias (deg/h)'); legend('X', 'Y', 'Z');