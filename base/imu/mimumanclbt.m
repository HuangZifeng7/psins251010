function [eb, Ka, db, imu1] = mimumanclbt(imu, t, g0, isfig)
% MIMU manual calibration.
%
% Prototype: [eb, Ka, db, imu1] = mimumanclbt(imu, t, g0, isfig)
% Inputs: imu - IMU before calibration = [wm, vm, t]
%         t - [static0,static1; x_up,x_down; y_up,y_down; z_up,z_down]
%         g0 - reference gravity
%         isfig - figure flag
% Outputs: eb, Ka, db - calibration results
%          imu - IMU output after calibration 
%
% See also  imuclbt, sysclbtMEMS, imuscale.

% Copyright(c) 2009-2022, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi'an, P.R.China
% 12/08/2022
global glv
    if nargin<4, isfig=1; end
    if nargin<3, g0=glv.g0; end
    ts = diff(imu(1:2,end));
    mm = zeros(7,7);
    for k=1:length(t)
        mm(k,:) = mean(datacut(imu,t(k,1),t(k,2)));
    end
    eb = mm(1,1:3)'/ts;
    xp = mm(2,4)/ts; xn = mm(3,4)/ts;
    yp = mm(4,5)/ts; yn = mm(5,5)/ts;
    zp = mm(6,6)/ts; zn = mm(7,6)/ts;
    Ka = diag([xp-xn; yp-yn; zp-zn]/2/g0);
    db = [xp+xn; yp+yn; zp+zn]/2;
    if nargout>3
        imu1 = [imu(:,1:3)-eb'*ts, imu(:,4:6)*Ka-db'*ts, imu(:,end)];
    end
    if isfig
        imuplot(imu1, imu);
    end