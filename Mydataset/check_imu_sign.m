% check_imu_sign.m
clear; clc; close all;
glvs;

% 1. 加载数据
try
    raw_imu = load('ICM20602.txt'); 
catch
    error('找不到 ICM20602.txt');
end

% IMU列: [Time, Gx, Gy, Gz, Ax, Ay, Az] (根据README)
% 取前 1000 个点（通常是静止或低速）
sample_data = raw_imu(1:1000, :);
ts = mean(diff(sample_data(:,1)));

% 计算平均加速度 (转为 m/s?)
% 原始数据是增量(m/s)，除以 ts 得到加速度
acc_mean = mean(sample_data(:, 5:7)) / ts;

fprintf('--- IMU Z轴符号检查 ---\n');
fprintf('Acc X mean: %.4f m/s^2\n', acc_mean(1));
fprintf('Acc Y mean: %.4f m/s^2\n', acc_mean(2));
fprintf('Acc Z mean: %.4f m/s^2\n', acc_mean(3));

if acc_mean(3) > 5
    fprintf('?? 警告：Z轴加速度为 正值 (+g)。\n');
    fprintf('   PSINS (FRD系) 通常期望静止时 Z轴为 负值 (-g)。\n');
    fprintf('   -> 建议在生成 mat 时给加速度取反，或者检查坐标系定义。\n');
elseif acc_mean(3) < -5
    fprintf('? 正常：Z轴加速度为 负值 (-g)。\n');
    fprintf('   如果此时还发散，说明初始姿态(Roll/Pitch)可能给错了。\n');
else
    fprintf('? 异常：Z轴数值不对 (既不是+g也不是-g)，请检查单位。\n');
end

figure; 
plot(sample_data(:,1), sample_data(:,7)/ts); 
title('Raw Acc Z (m/s^2)'); grid on;