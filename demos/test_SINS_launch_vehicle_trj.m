% Trajectory generation for launch vehicle simulation.
% See also  test_SINS_launch_vehicle_trj.
% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 27/04/2025
glvs
ts = 0.01;       % sampling interval
avp0 = [[85;0;45]*glv.deg; [0;0;0]; [29*glv.deg; 106*glv.deg; 450]]; % init avp
% trajectory segment setting
xxx = [];
seg = trjsegment(xxx, 'init',         0);
seg = trjsegment(seg, 'uniform',      10);
seg = trjsegment(seg, 'accelerate',   100, xxx, 1);
seg = trjsegment(seg, 'uniform',      10);
seg = trjsegment(seg, 'headdown',     30, 85/30);
seg = trjsegment(seg, 'uniform',      100);
% generate, save & plot
trj = trjsimu(avp0, seg.wat, ts, 1);
trjfile('trj10ms.mat', trj);
insplot(trj.avp);
imuplot(trj.imu);

