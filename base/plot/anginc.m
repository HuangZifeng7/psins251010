function ang = anginc(imu, t1, t2, xyz)
% Angular increment of IMU;
%
% Prototype: ang = anginc(imu, t1, t2, xyz)
% Inputs: imu - IMU array.
%         t1,t2 - begin/end time.
%         xyz - =-1/2/3 for x/y/z axis
% Output: ang - angular increment in rad
%
% See also  accmean, imumeanplot.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 01/10/2025
    if nargin<4, xyz=1:3; end
    ang = sum(datacut(imu(:,[xyz,end]),t1,t2))';
    ang = ang(xyz);
