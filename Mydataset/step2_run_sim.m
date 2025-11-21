% demo_Final_Victory.m
clear; clc; close all;
glvs;

%% 1. 准备数据 (既然数据是好的，我们只需确保单位是弧度)
% 加载数据
try
    load('ICM20602_Ready.mat');
catch
    error('请先运行 step1_make_mat.m 生成数据！');
end

% ? 二次确认：如果之前的 step1 里没转弧度，这里补救一下
% 检查一下 Latitude 大小，如果大于 1 (比如 30)，说明是度，要转弧度
if abs(gnss(1,1)) > 1.5 
    disp('检测到 GNSS 是度，正在转为弧度...');
    gnss(:,1:2) = gnss(:,1:2) * glv.deg;
    avp0(1:3) = avp0(1:3) * glv.deg; % 如果 avp0 也是度的话
end

%% 2. 设置参数 (稳健模式)
% 初始化
ins = insinit(avp0, ts);
ins.lever = ins_lever; 

% 误差参数 (MEMS 标准配置)
% 初始误差给大一点，让滤波器自己去收敛
avperr = avperrset([1;1;10], 1, 10); 
imuerr = imuerrset(0.5, 200, 0.5, 10000); 

% 关键：放宽观测噪声，防止炸飞
Rk = [1.0; 1.0; 3.0]; 
Pmin = [avperrset([0.2,1.0],0.01,0.2); gabias(0.1, [100,1000]); [0.01;0.01;0.01]; 0.001].^2;

%% 3. 解算 (关闭花里胡哨的功能)
disp('怪兽出击：开始解算...');

t_run = 600; % 跑600秒
max_k = min(length(imu), fix(t_run/ts));
imu_run = imu(1:max_k, :);

% 参数: 0=不画图, 'avp'=只修状态不修杆臂
% Acc不用取反，因为刚才检查了 Z 是 -9.7，是正确的！
[avp, xkpk, zkrk, sk, ins1, kf] = sinsgps(imu_run, gnss, ins, avperr, imuerr, ins.lever, ts, Rk, Pmin, 0, 'avp');

%% 4. 成果展示
if ~isempty(avp)
    disp('解算成功！正在绘图...');
    
    % 1. 轨迹图
    figure('Name', 'Victory Trajectory', 'Color', 'w');
    plot(avp(:,2)/glv.deg, avp(:,1)/glv.deg, 'b-', 'LineWidth', 2); hold on;
    
    % 提取 GNSS 对比
    t_sim = avp(:,end);
    valid_gnss = gnss(gnss(:,4)>=t_sim(1) & gnss(:,4)<=t_sim(end), :);
    plot(valid_gnss(:,2)/glv.deg, valid_gnss(:,1)/glv.deg, 'r--', 'LineWidth', 1.5);
    
    legend('GINS (Blue)', 'GNSS (Red)');
    title('Trajectory Comparison'); xlabel('Lon'); ylabel('Lat');
    grid on; axis equal;
    
    % 2. 看看高度是不是稳住了 (重点检查)
    figure('Name', 'Height Check', 'Color', 'w');
    plot(t_sim, avp(:,3)); hold on;
    plot(valid_gnss(:,4), valid_gnss(:,3), 'r--');
    legend('INS Height', 'GNSS Height');
    title('Height Check (Should be stable)'); grid on;
    
    % 3. 姿态
    insplot(avp);
else
    error('解算失败，AVP 为空');
end