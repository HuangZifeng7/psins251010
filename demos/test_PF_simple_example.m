% A simple example test for particle Kalman filter (PF).
% State space models:
%  参见：捷联惯导算法与组合导航原理（第二版）练习题-51（航空炸弹下落）
% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 12/02/2025
%% 轨迹参数设置
ts = 0.1;  T = 20;  len = T/ts;  dyfig = 0;   % 采用间隔，仿真时长，动态画粒子变化图标志
g = 9.8;  rho = 0.5;  % 重力，阻尼系数
L = 100;  H = 100;  % 雷达距离，炸弹初始高度
Qk = 0.2*ts;   sQ = chol(Qk);   w = sQ*randn(1,len);  % 系统噪声
Rk = diag([3,(1*pi/180)])^2;   sR = chol(Rk);   r = sR*randn(2,len);  % 量测噪声
X = zeros(2,len);  Z = zeros(2,len);  t = zeros(1,len);
X(:,1) = [H; 0];  Z(:,1) = [sqrt(L^2+X(1,1)^2); atan2(X(1,1),L)];  t(1) = ts;
%% 粒子滤波设置
N = 1000;  % 粒子数
Xk = X;  Pk = X;
dX0 = [10; 1];  Xk(:,1) = Xk(:,1)+dX0;  Pk(:,1) = Pk(:,1)+dX0.^2;  
particle = repmat(X(:,1),1,N) + diag(dX0)*randn(2,N);  % 初始化粒子
%% 粒子滤
for k = 2:len
    % 状态与量测仿真
    a = g - rho*X(2,k-1)^2;
    X(:,k) = X(:,k-1) + [-X(2,k-1)*ts-a*ts^2/2; a*ts+w(k)];  t(k)=k*ts;
    Z(:,k) = [sqrt(L^2+X(1,k)^2); atan2(X(1,k),L)] + r(:,k);
    % 粒子预测
    a = g - rho*particle(2,:).^2;
    particle = particle + [-particle(2,:)*ts-a*ts^2/2; a*ts+sQ*randn(1,N)];   
    Zpar = [sqrt(L^2+particle(1,:).^2); atan(particle(1,:)/L)] + sR*randn(2,N);
    if dyfig==1, figure(1); subplot(211), plot(particle(1,:)); subplot(212), plot(particle(2,:)); end
    % 权重计算
    zerr = [Z(1,k)-Zpar(1,:); Z(2,k)-Zpar(2,:)];
    Rzerr = Rk\zerr;  zRz=zerr(1,:).*Rzerr(1,:)+zerr(2,:).*Rzerr(2,:);
    weight = exp(-0.5*zRz);  weight = weight/sum(weight);
    % 重采样
    particle = particle(:,randsample(N,N,true,weight));
    % 状态估计及方差
    Xk(:,k) = mean(particle,2);  Pk(:,k) = diag(cov(particle'));
end
%% 作图
myfigure;
subplot(321), plot(t, Z(1,:)), grid on, xlabel('t / s'); ylabel('range / m');
subplot(322), plot(t, Z(2,:)*180/pi), grid on, xlabel('t / s'); ylabel('angle / \circ');
subplot(323), plot(t, [X(1,:); Xk(1,:)]), grid on, xlabel('t / s'); ylabel('H / m');
subplot(324), plot(t, [X(2,:); Xk(2,:)]), grid on, xlabel('t / s'); ylabel('V / m/s');
subplot(325), plot(t, [Xk(1,:)-X(1,:);sqrt(Pk(1,:))]), grid on, xlabel('t / s'); ylabel('Herr / m');
subplot(326), plot(t, [Xk(2,:)-X(2,:);sqrt(Pk(2,:))]), grid on, xlabel('t / s'); ylabel('Verr / m/s');
