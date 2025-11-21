% demo_Robust_Final_Rescue.m
% 专治：RCOND=NaN, 初始发散
% 核心修复：强制检查并转换 GPS 单位 (Deg -> Rad)
clear; clc; close all;
glvs;

%% 1. 数据加载与“防呆”预处理
ts = 0.01;
fprintf('1. 加载数据...\n');
try
    load mpu6500gsp.mat; 
catch
    error('找不到文件 mpu6500gsp.mat');
end

% 截取前 1500 秒
t_end = 1500;
idx_limit = min(length(imu), fix(t_end/ts));
imu = imu(1:idx_limit, :);
gps = gps(gps(:,end) <= imu(end,end), :);

% =======================================================
% ? 核心修复：单位自动检测与转换
% =======================================================
% 检查纬度是否 > 1.6 (pi/2)。如果是，说明是度，必须转弧度！
if abs(gps(1,1)) > 1.6
    fprintf('?? 检测到 GPS 单位是 度(Deg)，正在强制转为 弧度(Rad)...\n');
    gps(:,1:2) = gps(:,1:2) * glv.deg; % 转弧度
else
    fprintf('? 检测到 GPS 单位疑似 弧度(Rad)，保持不变。\n');
end
% =======================================================

% ? 注入干扰 (t=600~620s, 误差50m)
t_bad_start = 600;
t_bad_end = 620;
gps_poisoned = gps; 
bad_idx = gps(:,end) >= t_bad_start & gps(:,end) <= t_bad_end;
num_bad = sum(bad_idx);
if num_bad > 0
    % 纬度经度加误差 (注意现在已经是弧度了，要除以 Re)
    err = [50/glv.Re, 50/glv.Re, 50]; 
    gps_poisoned(bad_idx, 1:3) = gps_poisoned(bad_idx, 1:3) + repmat(err, num_bad, 1);
    fprintf('   已注入 50m 干扰 (t=600~620s)\n');
end

%% 2. 初始化 (使用更稳健的参数)
% 使用第一个 GPS 点初始化位置
avp0 = [[0;0;-100]*glv.deg; 0;0;0; gps(1, 1:3)']; 

% 这里的参数稍微给大一点，MPU6500 比较飘
avperr = avperrset([1;1;10], 1, 10); 
imuerr = imuerrset(0.5, 200, 0.5, 10000); 
% P阵下限
Pmin = [avperrset([0.2,1.0],0.01,0.2); gabias(0.1, [100,1000]); [0.01;0.01;0.01]; 0.001].^2;
% 观测噪声 R (位置精度 10m)
Rk_val = poserrset(10);

% 强制定义模型
psinstypedef(153);

%% 3. 双路解算 (完全相同的循环结构，只差一句抗差逻辑)
modes = {'Standard', 'Robust'};
res_avp = cell(1,2);

for m = 1:2
    mode = modes{m};
    fprintf('   -> 运行: %s ...\n', mode);
    
    % 重置 INS 和 KF
    ins = insinit(avp0, ts);
    ins.lever = [0;0;0]; % 暂时忽略杆臂求稳
    
    kf = kfinit(ins, avperr, imuerr, Rk_val);
    kf.Pmin = Pmin;
    kf.Rk = diag(Rk_val.^2);
    R_init = kf.Rk; % 备份初始 R
    
    n_state = size(kf.Pxk, 1); % 获取状态维度
    
    gps_ptr = 1;
    avp_hist = zeros(length(imu), 10);
    
    % 等待条
    bar = waitbar(0, [mode ' Running...']);
    
    for k = 1:length(imu)
        % INS
        ins = insupdate(ins, imu(k,:));
        % KF 预测
        kf.Phikk_1 = kffk(ins);
        kf = kfupdate(kf);
        
        % GNSS 更新
        curr_t = imu(k,end);
        if gps_ptr <= size(gps_poisoned, 1) && curr_t >= gps_poisoned(gps_ptr, end)
            
            pos_meas = gps_poisoned(gps_ptr, 1:3)';
            Z = pp2vn(ins.pos, pos_meas); % 残差
            
            % ==================================
            % ?? 抗差逻辑 (防止 NaN 的关键)
            % ==================================
            update_flag = 1;
            
            if strcmp(mode, 'Robust')
                % 阈值 25m
                innov = norm(Z);
                if innov > 25.0
                    % 策略：如果误差太大，直接不更新 (Gating)
                    % 这比膨胀 R 阵更安全，绝对不会导致矩阵病态
                    update_flag = 0; 
                end
            end
            % ==================================
            
            if update_flag
                kf.Hk = zeros(3, n_state);
                kf.Hk(1:3, 7:9) = eye(3);
                kf = kfupdate(kf, Z);
                % 反馈
                [kf, ins] = kffeedback(kf, ins, 1, 'avp');
            end
            
            gps_ptr = gps_ptr + 1;
        end
        
        avp_hist(k,:) = [ins.avp', curr_t];
        if mod(k, 5000)==0, waitbar(k/length(imu), bar); end
    end
    close(bar);
    res_avp{m} = avp_hist;
end

%% 4. 绘图 (下班图)
fprintf('3. 绘图...\n');
avp_std = res_avp{1};
avp_rob = res_avp{2};

% 恢复为度用于画图
pos_std = avp_std(:, 7:9); pos_std(:,1:2) = pos_std(:,1:2)/glv.deg;
pos_rob = avp_rob(:, 7:9); pos_rob(:,1:2) = pos_rob(:,1:2)/glv.deg;
bad_gps = gps_poisoned(bad_idx, :); bad_gps(:,1:2) = bad_gps(:,1:2)/glv.deg;

figure('Name', 'Final Fix', 'Color', 'w');

% 轨迹
subplot(2,1,1); hold on; grid on;
plot(bad_gps(:,2), bad_gps(:,1), 'rx', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(pos_std(:,2), pos_std(:,1), 'k--', 'LineWidth', 1);
plot(pos_rob(:,2), pos_rob(:,1), 'b-', 'LineWidth', 2);
legend('Interference (50m)', 'Standard KF', 'Robust KF');
title('Trajectory (Deg vs Rad Fixed)'); xlabel('Lon'); ylabel('Lat');

% 误差消除
subplot(2,1,2); hold on; grid on;
% 简单计算距离差
diff_m = sqrt(sum((avp_std(:,7:8) - avp_rob(:,7:8)).^2, 2)) * glv.Re;
plot(avp_rob(:,end), diff_m, 'g-', 'LineWidth', 1.5);
yl = ylim;
patch([t_bad_start t_bad_end t_bad_end t_bad_start], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
title('Optimization Benefit (m)'); xlabel('Time');

fprintf('完成！如果不出意外，这次轨迹应该正常了。\n');