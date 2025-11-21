function  [yxy, yxy0] = yxyrnd(ts, T, vel, n, isfig)
% Random trajectory yaw/X/Y generator, with 0-pitch/0-roll/0-Z and
%  constant velocity.

% Prototype: yxy = yxyrnd(ts, T, vel, n, isfig)
% Inputs: ts - sampling time interval
%         T - total time
%         vel - constant velocity
%         n - random way-point
%         isfig - figure flag
% Output: yxy = [yaw, x, y, t]
%
% Example
%   yxy = yxyrnd(0.01, 100, 10, 10);
%
% See also  imupolar_rnd, ap2avp, avptrans, avp2imu_ecef.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 14/12/2024
global glv
    if nargin<5, isfig = 1; end
    if nargin<4, n = 10; end;    if n<3, n=3; end
    if nargin<3, vel = 10; end
    if nargin<2, T = 100; end
    % to generate yaw
    if length(n)>1  % if user defined points n = [x,y]
        x = cumsum(n(:,1)); y = cumsum(n(:,2));
    else
        x = cumsum(randn(n,1));  y = cumsum(randn(n,1));
    end
    dis = sum(normv([x,y])); s = T*vel/dis;   x = x*s;  y = y*s;
    len = length(x);
    t = (0:len-1)'*T/(len-1);
    t1 = (0:ts:t(end)-ts)';
    x1 = spline(t,x,t1);  y1 = spline(t,y,t1);
    t2 = t1+ts/100;
    x2 = spline(t,x,t2);  y2 = spline(t,y,t2);
    yaw = atan2(y2-y1, x2-x1);
    % to generate x/y with constant velocity
    x = cumsum(-vel*ts*sin(yaw));  y = cumsum(vel*ts*cos(yaw));  % north point to y+
    x = x-mean(x); y = y-mean(y);
    if isfig==1
        myfig
        subplot(2,2,[1,3]), plot(x,y); xygo('x / m', 'y / m'); plot(x(1),y(1), 'or'); axis equal;
        subplot(2,2,2), plot(t1,yaw/glv.deg); xygo('y');
        subplot(2,2,4), plot(t1,[x,y]); xygo('xy / m');  title(sprintf('vel=%.2fm/2', vel));
    end
    yxy = [yaw, x, y, t1];
    yxy0 = [[interp1(t1,yaw,t), interp1(t1,x,t), interp1(t1,y,t), t]; yxy(end,:)];
