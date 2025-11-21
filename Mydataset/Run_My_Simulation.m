%% 针对 WHU-i2Nav GINS 数据集的专用转换脚本
% 对应 README 文档格式，严格转换单位和坐标系
% Author: Generated based on your provided README
clc; clear; close all;
glvs; % 加载 PSINS 全局变量

%% 1. 文件读取 (请根据实际文件名修改)
disp('正在读取数据...');
% 假设 IMU 文件名为 imu.txt, 参考真值文件名为 ref.nav
raw_imu = load('my_imu.txt');  
raw_nav = load('my_ref.nav'); 

disp('数据读取完成，开始格式转换...');

%% 2. 处理 IMU 数据
% README 定义: 
% Col 1: Time (s)
% Col 2-4: Angle Increment (rad) [FRD] -> 角度增量
% Col 5-7: Velocity Increment (m/s) [FRD] -> 速度增量

% 计算采样时间 ts
ts = raw_imu(2,1) - raw_imu(1,1); 
disp(['检测到 IMU 采样间隔 ts = ', num2str(ts), ' s']);

% PSINS 的 imu 矩阵通常需要 [Gyro(rad/s), Acc(m/s^2)] 
% 但对于增量数据，我们需要除以 ts 还原成“等效速率/比力”，或者在算法里指明是增量。
% 这里我们将其转换为标准的 Rate/Force 形式，这是最通用的做法。

imu = zeros(length(raw_imu), 7);

% [Gx, Gy, Gz] = 角度增量(rad) / ts -> rad/s
% 数据集已经是 FRD 系，PSINS 也是 FRD 系，【绝对不要】做 imurfu 旋转
imu(:, 1:3) = raw_imu(:, 2:4) / ts; 

% [Ax, Ay, Az] = 速度增量(m/s) / ts -> m/s^2
imu(:, 4:6) = raw_imu(:, 5:7) / ts; 

% 时间轴
imu(:, 7) = raw_imu(:, 1);

%% 3. 处理 Reference (Nav) 真值数据
% README 定义:
% Col 1: Week
% Col 2: Time (s)
% Col 3: Lat (deg), Col 4: Lon (deg), Col 5: Hgt (m)
% Col 6: Vn (m/s),  Col 7: Ve (m/s),  Col 8: Vd (m/s)
% Col 9: Roll(deg), Col 10: Pitch(deg), Col 11: Yaw(deg)

% PSINS 标准 avp 顺序: [Att(3), Vel(3), Pos(3)]
% 单位要求: Att(rad), Vel(m/s), Pos(Lat/Lon rad, Hgt m)

avp_ref = zeros(length(raw_nav), 10); % 第10列放时间

% --- 姿态 Att (Roll, Pitch, Yaw) ---
% 原始是 deg -> 转 rad
avp_ref(:, 1:3) = raw_nav(:, 9:11) * glv.deg;

% --- 速度 Vel (Vn, Ve, Vd) ---
% 原始是 m/s -> 保持不变
avp_ref(:, 4:6) = raw_nav(:, 6:8);

% --- 位置 Pos (Lat, Lon, Hgt) ---
% 原始 Lat/Lon 是 deg -> 转 rad
avp_ref(:, 7) = raw_nav(:, 3) * glv.deg; % Lat
avp_ref(:, 8) = raw_nav(:, 4) * glv.deg; % Lon
avp_ref(:, 9) = raw_nav(:, 5);           % Height

% --- 时间 ---
avp_ref(:, 10) = raw_nav(:, 2);

%% 4. 生成用于组合导航的 GPS 观测数据 (用于仿真输入)
% 提取 Nav 中的 Pos 和 Vel 作为模拟的 GPS 观测值
% 格式: [Vn, Ve, Vd, Lat, Lon, Hgt, Time]
gps_data = [avp_ref(:, 4:9), avp_ref(:, 10)];

%% 5. 保存为 .mat 文件
% 这里的变量命名非常关键，要配合 PSINS 的习惯
save('whu_gins_data.mat', 'imu', 'avp_ref', 'gps_data', 'ts');

disp('转换成功！文件已保存为 whu_gins_data.mat');
disp('包含变量: imu (FRD), avp_ref (真值), gps_data (观测), ts');

%% 6. 简单验证画图 (确保数据没飞)
figure;
subplot(2,1,1); 
plot(avp_ref(:,8)/glv.deg, avp_ref(:,7)/glv.deg); grid on;
title('真值轨迹 (Lon/Lat)'); xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
subplot(2,1,2);
plot(imu(:,7), imu(:,1:3)*glv.deg); grid on;
title('IMU 陀螺仪原始数据 (deg/s)'); xlabel('Time (s)');