% step1_make_mat.m
% 作用：清洗数据，统一为【弧度/米/秒】标准格式，剔除时间不对齐的数据
clear; clc; glvs;

disp('1. 读取原始文件...');
try
    raw_imu = load('ICM20602.txt'); 
    tr = load('truth.nav');
catch
    error('找不到文件！请确保 ICM20602.txt 和 truth.nav 在当前目录下');
end

% --- 1. 处理 IMU ---
% 原始列: [Time, dAng_x, dAng_y, dAng_z, dVel_x, dVel_y, dVel_z]
% PSINS标准: [Gx, Gy, Gz, Ax, Ay, Az, Time]
imu = raw_imu(:, [2:7, 1]); 

% 强制设定采样间隔 (避免因数据丢帧导致计算出的 dt 不准)
ts = 0.005; % 200Hz
disp(['   IMU 采样率设定为: ', num2str(1/ts), ' Hz']);

% --- 2. 处理 GNSS ---
% truth.nav 格式: [Week, Time, Lat(deg), Lon(deg), Hgt(m), ...]
gnss_full = tr(:, [2, 3, 4, 5]); 

% 截取时间：只保留 IMU 时间范围内的数据
t_start = imu(1,end);
t_end   = imu(end,end);
idx = gnss_full(:,1) > (t_start + 2.0) & gnss_full(:,1) < (t_end - 2.0);
gnss_data = gnss_full(idx, :);
gnss_data = gnss_data(1:100:end, :); % 降采样到 1Hz，避免相关性噪声

% ===> 核心：GNSS 必须转为 弧度 (Radians) <===
gnss = zeros(size(gnss_data,1), 4);
gnss(:,1) = gnss_data(:,2) * glv.deg; % Lat -> Rad
gnss(:,2) = gnss_data(:,3) * glv.deg; % Lon -> Rad
gnss(:,3) = gnss_data(:,4);           % Hgt -> m
gnss(:,4) = gnss_data(:,1);           % Time

% 添加 10cm 级别的白噪声模拟 RTK
gnss(:,1) = gnss(:,1) + (0.1/6378137)*randn(size(gnss,1),1);
gnss(:,2) = gnss(:,2) + (0.1/6378137)*randn(size(gnss,1),1);
gnss(:,3) = gnss(:,3) + 0.2*randn(size(gnss,1),1);

% --- 3. 初始对准参数 ---
% 找到 IMU 开始时刻最近的真值
[~, idx_init] = min(abs(tr(:,2) - t_start));
init_tr = tr(idx_init, :);

% 构造初始 AVP (全部转弧度)
% Truth: 9=Roll, 10=Pitch, 11=Yaw (deg)
att0 = [init_tr(10); init_tr(9); init_tr(11)] * glv.deg; 
vel0 = init_tr(6:8)'; % Vn, Ve, Vd
pos0 = [init_tr(3)*glv.deg; init_tr(4)*glv.deg; init_tr(5)]; % Lat, Lon, Hgt
avp0 = [att0; vel0; pos0];

% 固定杆臂 (来自 README)
ins_lever = [-0.073; 0.302; 0.087]; 

save('ICM20602_Ready.mat', 'imu', 'gnss', 'avp0', 'ts', 'ins_lever');
disp('数据准备完成：ICM20602_Ready.mat');