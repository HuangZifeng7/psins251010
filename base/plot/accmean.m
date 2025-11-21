function acc = accmean(imu, t1, t2, xyz)
% Acc mean specific force of IMU;
%
% Prototype: acc = accmean(imu, t1, t2, xyz)
% Inputs: imu - IMU array.
%         t1,t2 - begin/end time.
%         xyz - =-1/2/3 for x/y/z axis
% Output: acc - mean specific force in m/s^2
%
% See also  anginc, imumeanplot.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 01/10/2025
    [m,n] = size(imu);
    if n<6, n=0; else n=3; end
    if nargin<4, xyz=1:3; end
    acc = mean(datacut(imu(:,[xyz+n,end]),t1,t2))';
    acc = acc(xyz)/diff(imu(1:2,end));
