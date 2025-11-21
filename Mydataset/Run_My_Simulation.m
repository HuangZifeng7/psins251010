% ...existing code...
disp('开始读取数据...');

% --- 1) 检查文件是否存在 ---
imu_fn = 'my_imu.txt';
nav_fn = 'my_ref.nav';
if exist(imu_fn,'file')~=2, error('找不到 %s', imu_fn); end
if exist(nav_fn,'file')~=2, error('找不到 %s', nav_fn); end

% --- 2) 读取原始数据（更通用的 readmatrix 容错） ---
try
    raw_imu = load(imu_fn);   % 兼容纯数值txt
catch
    raw_imu = readmatrix(imu_fn);
end
try
    raw_nav = load(nav_fn);
catch
    raw_nav = readmatrix(nav_fn);
end

% --- 3) 查找 IMU 时间列（单调递增的列） ---
[nr_imu, nc_imu] = size(raw_imu);
timecol_imu = find(arrayfun(@(c) all(diff(raw_imu(:,c))>0), 1:nc_imu), 1);
if isempty(timecol_imu)
    error('无法自动识别 IMU 时间列：raw_imu 时间不单调。请检查文件格式。');
end
imu_time = raw_imu(:, timecol_imu);

% --- 4) 根据列数和常见约定挑 gyro/acc 列 ---
% 最常见格式： [time, gx, gy, gz, ax, ay, az, ...]
cols = setdiff(1:nc_imu, timecol_imu);
if length(cols) >= 6
    gyr_cols = cols(1:3);
    acc_cols = cols(4:6);
else
    error('IMU 列数不足，期望至少包含 time,gx,gy,gz,ax,ay,az');
end
imu_gyr = raw_imu(:, gyr_cols);  % 可能是增量或速率
imu_acc = raw_imu(:, acc_cols);

% --- 5) 估计采样间隔 ts 并做简单检验 ---
dt = diff(imu_time);
if any(dt<=0), error('IMU 时间不单调或含重复'); end
ts = median(dt);
disp(['检测到 IMU 采样间隔 ts = ', num2str(ts), ' s (median)']);

% --- 6) 判定 imu_type (heuristic) 并保证保存为 delta(增量) 格式 ---
% 经验阈值：若角数据平均绝对值 < 1e-1 则很有可能是 delta-angle（rad）
gyr_mean = mean(abs(imu_gyr(:)));
acc_mean = mean(abs(imu_acc(:)));
if gyr_mean < 1e-1 && acc_mean < 1e-1
    imu_type = 'delta';    % 已经是每步增量（不乘 dt）
    imu = [imu_gyr, imu_acc, imu_time];
else
    imu_type = 'rate';     % 假定是角速率 / 加速度 -> 转增量
    imu = [imu_gyr*ts, imu_acc*ts, imu_time];
    disp('注意：检测到 imu 看起来像 rate，脚本已把其乘以 ts 转为 delta 格式');
end
disp(['imu_type 判定为: ', imu_type, ' (已保存为 delta 格式)']);
disp(['imu 大小：', num2str(size(imu,1)), ' x ', num2str(size(imu,2))]);

% --- 7) 解析 raw_nav -> avp_ref / gps_data ---
[nr_nav, nc_nav] = size(raw_nav);
% 自动识别时间列（单调递增）
timecol_nav = find(arrayfun(@(c) all(diff(raw_nav(:,c))>0), 1:nc_nav), 1);
if isempty(timecol_nav)
    % 若找不到单调列，尝试第一列为时间
    timecol_nav = 1;
end
% 尝试识别 lat/lon/height 列：值范围判断
latcol = find(all(raw_nav >= -90 & raw_nav <= 90, 1), 1); % deg 范围
loncol = find(all(raw_nav >= -180 & raw_nav <= 180, 1) & (1:nc_nav)~=latcol, 1);
hcol   = find((1:nc_nav)~=timecol_nav & (1:nc_nav)~=latcol & (1:nc_nav)~=loncol, 1);

% 若推断失败，退回默认顺序
if isempty(latcol) || isempty(loncol) || isempty(hcol)
    % 常见格式: [time lat lon h vn ve vd] 或 [lat lon h time vn ve vd]
    if nc_nav >= 7
        if timecol_nav == 1
            tcol = 1; latcol = 2; loncol = 3; hcol = 4; vncol = 5; vecol = 6; vdcol = 7;
        else
            tcol = 4; latcol = 1; loncol = 2; hcol = 3; vncol = 5; vecol = 6; vdcol = 7;
        end
    elseif nc_nav == 4
        % [time lat lon h]
        tcol = 1; latcol = 2; loncol = 3; hcol = 4;
        vncol = []; vecol = []; vdcol = [];
    else
        error('raw_nav 列数/格式无法自动识别，请手动检查文件 head');
    end
else
    % 确定时间列与速率列（简单分配）
    tcol = timecol_nav;
    % choose other columns by elimination (best-effort)
    % assume vn/ve/vd 最常在末尾
    other = setdiff(1:nc_nav, [tcol, latcol, loncol, hcol]);
    if length(other) >= 3
        vncol = other(1); vecol = other(2); vdcol = other(3);
    else
        vncol = []; vecol = []; vdcol = [];
    end
end

time_nav = raw_nav(:, tcol);
lat_deg = raw_nav(:, latcol);
lon_deg = raw_nav(:, loncol);
h_m  = raw_nav(:, hcol);
if ~isempty(vncol)
    vn = raw_nav(:, vncol); ve = raw_nav(:, vecol); vd = raw_nav(:, vdcol);
else
    vn = zeros(nr_nav,1); ve = zeros(nr_nav,1); vd = zeros(nr_nav,1);
end

% 转换为弧度/米并组装
lat_rad = lat_deg * pi/180;
lon_rad = lon_deg * pi/180;
% --- Robustness check: sometimes files have lat/lon order swapped or unexpected units
% If mean latitude magnitude is very small (e.g. < 1 deg) while longitude is large (e.g. > 5 deg)
% then assume lat/lon might be swapped in raw data and swap them.
if mean(abs(lat_deg)) < 1 && mean(abs(lon_deg)) > 5
    % likely swapped; swap lat/lon
    warning('检测到 raw_nav 中 lat/lon 似乎颠倒（lat 平均值很小但 lon 很大），已自动交换列。');
    tmp = lat_deg; lat_deg = lon_deg; lon_deg = tmp;
    lat_rad = lat_deg * pi/180; lon_rad = lon_deg * pi/180;
end
% 解析真实姿态（如果 raw_nav 提供 roll/pitch/yaw）
att0 = zeros(nr_nav,3);
% 许多数据集在后边会包含 roll/pitch/yaw（列索引可能是 9..11 或末尾几列）
if nc_nav >= 11
    % 常见顺序： [week, sow, lat, lon, h, vn, ve, vd, roll, pitch, yaw]
    % 先尝试按 9:11 为 roll/pitch/yaw
    rollcol = 9; pitchcol = 10; yawcol = 11;
    if all(abs(raw_nav(:, rollcol)) <= 360) && all(abs(raw_nav(:, pitchcol)) <= 360) && all(abs(raw_nav(:, yawcol)) <= 360)
        att0 = raw_nav(:, [rollcol, pitchcol, yawcol]) * pi/180; % deg->rad
    end
end
vel = [vn, ve, vd];
pos = [lat_rad, lon_rad, h_m];
avp_ref = [att0, vel, pos, time_nav];

% GNSS 格式： [Vn, Ve, Vd, Lat(rad), Lon(rad), Hgt(m), Time]
gps_data = [vn, ve, vd, lat_rad, lon_rad, h_m, time_nav];

% --- 8) 保存（带 timestamp，避免覆盖） ---
outfn = ['whu_gins_data_', datestr(now,'yyyymmdd_HHMMSS'), '.mat'];
save(outfn, 'imu', 'imu_type', 'avp_ref', 'gps_data', 'ts');
% 也写出固定文件名以供解算脚本直接加载
save('whu_gins_data.mat', 'imu', 'imu_type', 'avp_ref', 'gps_data', 'ts');
disp(['转换成功并保存为: ', outfn, ' (同时覆盖 whu_gins_data.mat)']);

% --- 9) 简单验证显示 ---
disp('变量信息：');
whos imu imu_type avp_ref gps_data ts
fprintf('imu_time 范围: [%g, %g]\n', imu(1,end), imu(end,end));
fprintf('gps_time 范围: [%g, %g]\n', gps_data(1,end), gps_data(end,end));
% 显示第一行样例，便于确认列映射是否正确
disp('IMU 首行示例:'); disp(imu(1,:));
disp('GPS 首行示例:'); disp(gps_data(1,:));

% ...existing code...
% (已在上面保存并生成 whu_gins_data.mat)
% ...existing code...

disp('转换成功！已生成并覆盖：whu_gins_data.mat');
disp('包含变量: imu (FRD delta), imu_type, avp_ref (真值), gps_data (观测), ts');

%% 6. 简单验证画图 (确保数据没飞)
figure;
subplot(2,1,1); 
plot(avp_ref(:,8)/glv.deg, avp_ref(:,7)/glv.deg); grid on;
title('真值轨迹 (Lon/Lat)'); xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
subplot(2,1,2);
plot(imu(:,7), imu(:,1:3)*glv.deg); grid on;
title('IMU 陀螺仪原始数据 (deg/s)'); xlabel('Time (s)');