% ========================================================================
% SINS/GPS 组合导航仿真：绝地反击版 (修正数值爆炸 + 修正索引错误)
% ========================================================================
clear; clc; close all;
glvs;
psinstypedef(153); 

%% 1. 数据加载
fprintf('1. 加载数据...\n');
try
    dd = load('data2.txt');
catch
    error('找不到 data2.txt，请确认文件路径！');
end

ts = 0.01; 
tt = dd(:,1);
imu = [[dd(:,2:4)*glv.dps, dd(:,5:7)*glv.g0]*ts, tt]; 

% --- 修正点 1: 明确 trj_avp 的列结构 ---
% Col 1-3: Att (姿态)
% Col 4-6: Vel (速度)
% Col 7-9: Pos (位置: Lat, Lon, Hgt)
% Col 10:  Time (时间)
trj_avp = [dd(:,8:10)*glv.deg, dd(:,11:13), dd(:,14:15)*glv.deg, dd(:,16), tt];
trj_avp(:,3) = yawcvt(trj_avp(:,3),'c360cc180');

% 原始 GPS 
gps_clean = [dd(:,17:19), dd(:,20:22), tt];
gps_clean = gps_clean(gps_clean(:,4)>0, :);

%% 2. 构造野值 (干扰环境)
gps_bad = gps_clean;
mid_idx = fix(length(gps_bad)/2);
bad_range = mid_idx : mid_idx+30; 
gps_bad(bad_range, 1:3) = gps_bad(bad_range, 1:3) + 30; % 注入30米偏差
fprintf('   场景：中间注入 30米 GPS 干扰。\n');

%% 3. 参数设置
[nn, ts, nts] = nnts(2, ts);
imuerr = imuerrset(0.03, 100, 0.001, 5);
davp0 = avperrset([1;1;10], 0.1, [1;1;3]); % 初始误差设小一点，保命要紧
rk = poserrset(10); 

%% 4. 运行算法 1: 标准 KF
fprintf('2. 运行标准 KF ...\n');
ins = insinit(avpadderr(trj_avp(1,1:9)', davp0), ts); 
kf = kfinit(ins, davp0, imuerr, rk);
kf.Pmin = [avperrset(0.01,1e-4,0.1); gabias(1e-3, [1,10])].^2;  
kf.pconstrain = 1; 

len = length(imu);
[avp_kf, xkpk_kf] = prealloc(fix(len/nn), 10, 2*kf.n+1);
ki = 1;
gps_ptr = 1;

for k=1:nn:len-nn+1
    k1 = k+nn-1;
    wvm = imu(k:k1, 1:6); t = imu(k1, end);
    ins = insupdate(ins, wvm); 
    
    kf.Phikk_1 = kffk(ins);    
    kf = kfupdate(kf);         
    
    if gps_ptr <= length(gps_bad) && abs(t - gps_bad(gps_ptr,end)) < nts
        posGPS = gps_bad(gps_ptr, 1:3)'; 
        kf = kfupdate(kf, ins.pos - posGPS, 'M'); 
        [kf, ins] = kffeedback(kf, ins, 1, 'avp');
        gps_ptr = gps_ptr + 1;
    end
    avp_kf(ki,:) = [ins.avp', t];
    ki = ki+1;
end
avp_kf(ki:end,:) = []; 

%% 5. 运行算法 2: 抗差 EKF (数值安全版)
fprintf('3. 运行抗差 EKF ...\n');
ins = insinit(avpadderr(trj_avp(1,1:9)', davp0), ts); 
kf = kfinit(ins, davp0, imuerr, rk);
kf.Pmin = [avperrset(0.01,1e-4,0.1); gabias(1e-3, [1,10])].^2;  
kf.pconstrain = 1;

[avp_rob, xkpk_rob] = prealloc(fix(len/nn), 10, 2*kf.n+1);
ki = 1;
gps_ptr = 1;
IGG_Threshold = 5.0; 

for k=1:nn:len-nn+1
    k1 = k+nn-1;
    wvm = imu(k:k1, 1:6); t = imu(k1, end);
    
    ins = insupdate(ins, wvm);
    
    kf.Phikk_1 = kffk(ins); 
    kf.px = ins;            
    kf = ekf(kf);           
    
    if gps_ptr <= length(gps_bad) && abs(t - gps_bad(gps_ptr,end)) < nts
        posGPS = gps_bad(gps_ptr, 1:3)';
        zk = ins.pos - posGPS; 
        
        % 1. 构造 H
        H = zeros(3, kf.n);      
        H(:, 7:9) = eye(3);      
        
        % 2. 计算马氏距离 (增加NaN检查)
        % 强制对称化 P 阵，防止数值不稳定
        kf.Pxk = (kf.Pxk + kf.Pxk')/2;
        P = kf.Pxk;
        R = kf.Rk; 
        S = H*P*H' + R;
        
        % 使用伪逆 pinv 防止矩阵奇异导致的崩溃
        try
            beta = zk' * (pinv(S) * zk); 
        catch
            beta = 0;
        end
        
        % 如果 beta 计算出来是 NaN，强制设为大数(视为野值)
        if isnan(beta), beta = 9999; end
        
        % 3. 抗差策略：检测到野值，直接跳过更新！(最安全的做法)
        if sqrt(beta) > IGG_Threshold
            % 策略：Do Nothing (相当于 R = 无穷大)
            % fprintf('跳过野值: t=%.2f\n', t);
        else
            % 只有正常点才更新
            kf = kfupdate(kf, zk, 'M');
        end
        
        [kf, ins] = kffeedback(kf, ins, 1, 'avp');
        gps_ptr = gps_ptr + 1;
    end
    
    avp_rob(ki,:) = [ins.avp', t];
    ki = ki+1;
end
avp_rob(ki:end,:) = [];

%% 6. 绘图 (修正索引错误)
fprintf('4. 正在绘图...\n');

% 检查数据
if isempty(avp_kf) || isempty(avp_rob)
    error('仿真结果为空，请检查数据输入。');
end

% 计算误差
% 插值对齐 (trj_avp 只有10列，时间在第10列)
t_ref = trj_avp(:, 10); 
t_kf  = avp_kf(:, end);
t_rob = avp_rob(:, end);

% 提取真值的位置 (Col 7,8,9)
ref_pos_kf = zeros(length(t_kf), 3);
ref_pos_rob = zeros(length(t_rob), 3);

for i=1:3
    ref_pos_kf(:,i)  = interp1(t_ref, trj_avp(:,6+i), t_kf, 'linear', 'extrap');
    ref_pos_rob(:,i) = interp1(t_ref, trj_avp(:,6+i), t_rob, 'linear', 'extrap');
end

% KF 位置在 7,8,9 列
err_kf_pos  = avp_kf(:,7:9)  - ref_pos_kf;
err_rob_pos = avp_rob(:,7:9) - ref_pos_rob;

% === 绘制误差图 ===
figure('Name', '位置误差对比', 'Color', 'w', 'NumberTitle', 'off');
titles = {'北向误差 (m)', '东向误差 (m)', '地向误差 (m)'};
for i=1:3
    subplot(3,1,i);
    plot(t_kf, err_kf_pos(:,i), 'b'); hold on;
    plot(t_rob, err_rob_pos(:,i), 'r', 'LineWidth', 1.5);
    ylabel(titles{i}); grid on;
    if i==1
        title('SINS/GPS 位置误差对比 (含野值)');
        legend('标准 KF', '抗差 EKF');
    end
    % 自动缩放 Y 轴，但排除极端的 NaN 或 Inf
    valid_idx = abs(err_kf_pos(:,i)) < 1000; 
    if any(valid_idx)
        ylim([-50, 50]); % 强制锁定在合理范围，发散的线让它飞出画面
    end
end
xlabel('时间 (s)');

% === 修正后的轨迹图 ===
figure('Name', '平面轨迹对比', 'Color', 'w');
% --- 修正点 2: trj_avp 的经纬度在第 8 和 7 列 (Lon, Lat) ---
% 之前你代码里写的是 12 和 11，那是错的
plot(trj_avp(:,8), trj_avp(:,7), 'k--', 'LineWidth', 1.5); hold on; 
plot(avp_kf(:,8), avp_kf(:,7), 'b'); 
plot(avp_rob(:,8), avp_rob(:,7), 'r');

legend('真值', '标准 KF', '抗差 EKF');
xlabel('经度'); ylabel('纬度'); grid on;
axis equal;

fprintf('大功告成！\n');