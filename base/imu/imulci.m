function imu = imulci(imu, typ)
% Trans R/F/U IMU direction to LCI F/U/R direction.
%
% Prototype: imu = imulci(imu, typ)
% Inputs: imu - IMU in.
%         typ - 1 for R/F/U to F/U/R, 0 for F/U/R to R/F/U
% Output: imu - IMU out.
%
% See also  imurfu.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 18/05/2025
    if nargin<2, typ=1; end
    if typ==1
        imu(:,1:6) = imu(:,[2,3,1, 5,6,4]);
    else
        imu(:,1:6) = imu(:,[3,1,2, 6,4,5]);
    end