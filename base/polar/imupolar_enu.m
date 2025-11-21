function [imu, avp0, avp] = imupolar_enu(pos0, vel0, ts, T, imuerr)
% SIMU sensor incremental outputs on polar filght.
%
% Prototype: [imu, avp0, avp] = imupolar_enu(pos0, vel0, ts, T, imuerr)
% Inputs: pos0 - initial pos0=[lat,lon,hgt]
%         vel0 - initial velocity forward, with header pointing to east(yaw=-90deg)
%         ts - SIMU sampling interval
%         T - total sampling simulation time
%         imuerr - SIMU error setting structure array from imuerrset
% Outputs: imu - gyro & acc incremental outputs
%          avp0 - initial avp in ENU-frame
%          avp - avp series in ENU-frame
%
% Example:
%   [imu, avp0, avp] = imupolar_enu(posset(89.5,0,0), 100, 1, 3000);
%   imuplot(imu);  insplot_polar(avp,'g');
%
% See also  imupolar_ecef, imupolar_grid, avp2imu_ecef, imustatic.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
global glv
    [RMh,clRNh] = RMRN(pos0);
    t = (0:ts:T)';  len=length(t);
    lat = repmat(pos0(1), len,1);
    lon = pos0(2)+t*vel0/clRNh;
    hgt = repmat(pos0(3), len,1);
    avp = avptrans([zeros(len,2), -pi/2*ones(len,1), zeros(len,3), lat,lon,hgt], 'n2e');  % grid-frame to ENU-frame, with grid yaw==0
    avp = ap2avp_ecef([avp,t]);
    [imu, avp0] = avp2imu_ecef(avp);
    avp = avptrans(avp, 'e2n');
    if exist('imuerr', 'var')
        imu = imuadderr(imu, imuerr);
    end
