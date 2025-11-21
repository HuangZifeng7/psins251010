%% WHU-i2Nav GINS 数据集 SINS/GNSS 组合导航解算脚本
% 适配数据: whu_gins_data.mat (由之前的脚本转换得到)
% 核心逻辑: 15维误差状态卡尔曼滤波 (松耦合)
% Time: 2025
clc; clear; close all;
glvs; 
psinstypedef(153); % 使用15维状态模型

%% 1. 加载数据
disp('正在加载数据...');
if ~exist('whu_gins_data.mat', 'file')
    error('找不到 whu_gins_data.mat，请先运行数据转换脚本！');
end
load('whu_gins_data.mat'); 
% 变量: imu(FRD, rad/s, m/s^2), avp_ref(真值), gps_data(观测), ts

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
    if gps_ptr <= size(gps_data, 1) && t_now >= gps_data(gps_ptr, end)
        
         % === 【新增探测代码】 ===
       % disp(['正在执行第 ', num2str(gps_ptr), ' 次 GPS 更新...']);
        %disp(['当前位置差 (INS-GPS): ', num2str(ins.pos' - gps_data(gps_ptr, 4:6))]);
        
        % 获取当前时刻的 GPS 观测值
        pos_gps = gps_data(gps_ptr, 4:6)'; % Lat, Lon, Hgt
        vn_gps  = gps_data(gps_ptr, 1:3)'; % Vn, Ve, Vd
        
        % 构造位置量测差值 Z = INS_pos - GPS_pos
        % PSINS 内部会自动处理杆臂效应 (ins.ant)
        Z = ins.pos - pos_gps;
        
        % 执行量测更新 (这里只用了位置更新，如果你想加速度更新也可以)
        kf = kfupdate(kf, Z, 'M');
        
        % [闭环反馈] 将估计的误差修正到 INS 状态中
        [kf, ins] = kffeedback(kf, ins, 1, 'avp');
        
        % 指针后移
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
avp_ref_interp = interp1(avp_ref(:,end), avp_ref(:,1:9), res_avp(:,end), 'linear');

% 调用 PSINS 绘图工具
% 1. 轨迹对比
figure;
plot(gps_data(:,5)/glv.deg, gps_data(:,4)/glv.deg, 'k.'); hold on;
plot(res_avp(:,5)/glv.deg, res_avp(:,4)/glv.deg, 'r'); 
legend('GNSS Observed', 'SINS/GNSS Integrated');
xlabel('Longitude (deg)'); ylabel('Latitude (deg)'); grid on;
title('2D Trajectory');

% 2. 误差曲线
avp_err = res_avp(:, 1:9) - avp_ref_interp;
% 姿态误差转为角分/度
avp_err(:,1:3) = avp_err(:,1:3) / glv.min; % 转换为角分 (arcmin)
avp_err(:,4:9) = avp_err(:,4:9);           % 速度位置保持 m/s, m

figure;
subplot(3,1,1); plot(res_avp(:,end), avp_err(:,1:3)); grid on;
title('Attitude Error (arcmin)'); legend('\phi', '\theta', '\psi');
subplot(3,1,2); plot(res_avp(:,end), avp_err(:,4:6)); grid on;
title('Velocity Error (m/s)'); legend('Vn', 'Ve', 'Vd');
subplot(3,1,3); plot(res_avp(:,end), avp_err(:,7:9)); grid on;
title('Position Error (m)'); legend('dN', 'dE', 'dD');

% 3. 陀螺仪零偏估计曲线
figure;
plot(res_xk(:,end), res_xk(:,10:12)/glv.dph); grid on;
title('Estimated Gyro Bias (deg/h)'); legend('X', 'Y', 'Z');