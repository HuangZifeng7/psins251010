function [imu, avp0, avp] = imupolar_grid(pos0, vel0, ts, T, imuerr)
% SIMU sensor incremental outputs on polar filght.
%
% Prototype: [imu, avp0, avp] = imupolar_grid(pos0, vel0, ts, T, imuerr)
% Inputs: pos0 - initial pos0=[lat,lon,hgt]
%         vel0 - initial velocity forward, with header pointing to polar
%         ts - SIMU sampling interval
%         T - total sampling simulation time
%         imuerr - SIMU error setting structure array from imuerrset
% Outputs: imu - gyro & acc incremental outputs
%          avp0 - initial avp in ENU-frame
%          avp - avp series in grid-frame
%
% Example:
%   [imu, avp0, avp] = imupolar_grid(posset(89.0,10,100), 100, 1, 2000);
%   imuplot_polar(imu);  insplot_polar(avp,'g');
%
% See also  imupolar_rnd, avp2imu_ecef, imupolar_ecef, imupolar_enu, imustatic, avpcmpplot_polar.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 02/12/2024
    RMh = RMRN(pos0);  xyz0 = blh2xyz(pos0);
    t = (0:ts:T)';  len=length(t);
    lat = pos0(1)+t*vel0/RMh;
    idx=lat>pi/2; lat(idx)=pi-lat(idx);
    lon = repmat(pos0(2), len,1);   lon(idx) = pos0(2)+pi;
    hgt = repmat(pos0(3), len,1);
    llh = [lat,lon,hgt];
    for k=1:10
        xyz = blh2xyzBatch(llh);
        if max(abs(xyz(:,2)-xyz0(2)))<0.001, break; end
        xyz(:,2) = xyz0(2); llh = xyz2llhBatch(xyz); llh(:,3) = pos0(3);
    end
    avp = avptrans([zeros(len,6), xyz], 'g2e');  % grid-frame to ECEF-frame, with grid yaw==0
    avp = ap2avp_ecef([avp,t]);
    [imu, avp0] = avp2imu_ecef(avp);
    avp = avptrans(avp, 'e2g');
    if exist('imuerr', 'var')
        imu = imuadderr(imu, imuerr);
    end

