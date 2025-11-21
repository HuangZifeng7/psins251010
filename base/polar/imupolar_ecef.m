function [imu, avp0, avp] = imupolar_ecef(pos0, vel0, ts, T, imuerr)
% SIMU sensor incremental outputs on polar filght.
%
% Prototype: [imu, avp0, avp] = imupolar_ecef(pos0, vel0, ts, T, imuerr)
% Inputs: pos0 - initial pos0=[lat,lon,hgt]
%         vel0 - initial velocity forward, with header pointing to polar
%         ts - SIMU sampling interval
%         T - total sampling simulation time
%         imuerr - SIMU error setting structure array from imuerrset
% Outputs: imu - gyro & acc incremental outputs
%          avp0 - initial avp in ENU-frame
%          avp - avp series in ECEF-frame
%
% Example:
%   [imu, avp0, avp] = imupolar_ecef(posset(89,30,0), 100, 1, 3000);
%   imuplot(imu);  insplot_polar(avp,'e');
%
% See also  imupolar_grid, imupolar_enu, avp2imu_ecef, imustatic.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
global glv
    RMh = RMRN(pos0);  ns=sign(pos0(1));  % north+ or south-
    t = (0:ts:T)';  len=length(t);
    lat = pos0(1)+t*(ns*vel0/RMh);
    if ns>0, idx=lat>pi/2; lat(idx)=pi-lat(idx);
    else,    idx=lat<-pi/2; lat(idx)=-pi-lat(idx); end
    lon = repmat(pos0(2), len,1);   lon(idx) = pos0(2)+pi;
    hgt = repmat(pos0(3), len,1);
    if ns==1, yaw=zeros(len,1); else, yaw=repmat(pi,len,1); end;  yaw(idx) = yaw(1)-pi;
    avp = avptrans([zeros(len,2), yaw, zeros(len,3), lat,lon,hgt],'n2e');  % ENU-frame to ECEF-frame
    avp = ap2avp_ecef([avp(:,[1:3,7:9]),t]);
    [imu, avp0] = avp2imu_ecef(avp);
    if exist('imuerr', 'var')
        imu = imuadderr(imu, imuerr);
    end

