function imuview(imu)
% View SIMU data in the Simulation Data Inspector.
%
% Prototype: imuview(imu)
% Inputs: imu - SIMU data, the last column is time tag
%          
% See also  imuplot, dataview.
%
% Example
%   imuview(imu);

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 03/09/2025
global glv
    ts = diff(imu(1:2,end));
    dps = glv.dps; g0 = glv.g0;
    Simulink.sdi.view; Simulink.sdi.clear;
    viewRun = Simulink.sdi.Run.create;
    for k=1:3
        tsk = timeseries(imu(:,k)/dps/ts, imu(:,end), 'Name', ['Gyro-', char('x'+k-1)]);
        add(viewRun,'vars',tsk);
    end
    for k=4:6
        tsk = timeseries(imu(:,k)/ts, imu(:,end), 'Name', ['Acc-', char('x'+k-4)]);
        add(viewRun,'vars',tsk);
    end
    % Simulink.sdi.close