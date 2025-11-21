function  [imu, avp0, avp, yxy0] = imupolar_rnd(ts, T, vel, hgt, L0, isfig)
% Polar random trajectory simulation, with 0-pitch/0-roll/0-Z and
%  constant velocity & height.

% Prototype: [imu, avp0, avp] = imupolar_rnd(ts, T, vel, hgt, L0, isfig)
% Inputs: ts - sampling time interval
%         T - total time
%         vel,hgt - constant velocity & height
%         L0 - mean latitude in (x+,0), default (0,0)
%         isfig - figure flag
% Outputs: imu - gyro & acc incremental outputs
%          avp0 - initial avp in ENU-frame
%          avp - avp series in grid-frame
%
% Example
%   [imu, avp0, avp] = imupolar_rnd(0.01, 100, 10, 200, 89.9999*glv.deg);  % imuplot(imu);
%   avpg = inspure_grid(imu, avp0);  avpcmpplot_polar(avp, avpg);
%   avpn = inspure_enu(imu, avp0);   avpcmpplot_polar(avp, avptrans(avpn,'n2g'));
%
% See also  imupolar_grid, yxyrnd, ap2avp, avptrans, avp2imu_ecef.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 14/12/2024
global glv
    if nargin<6, isfig=1; end
    if nargin<5, L0=90*glv.deg; end
    if nargin<4, hgt=0; end
    if nargin<3, vel=10; end
    if nargin<2, T=100; end
    if nargin<1, ts=0.01; end
    [yxy, yxy0] = yxyrnd(ts, T, vel, 10, 0);  [m,n] = size(yxy);
    yxy(:,2) = yxy(:,2)+(90*glv.deg-L0)*glv.Rp;
    yxy0(:,2) = yxy0(:,2)+(90*glv.deg-L0)*glv.Rp;
    avp = [zeros(m,2),yxy(:,1)-pi/2, zeros(m,3), yxy(:,2:3),repmat(glv.Rp,m,1),yxy(:,4)];  % -pi/2: grid-north point to x-
    llh=xyz2llhBatch(avp(:,7:9));  llh(:,3)=hgt;  avp(:,7:9)=blh2xyzBatch(llh);  % insplot_polar(avp, 'g');
    avpe = avptrans(avp, 'g2e');  avpe = ap2avp_ecef(avpe(:,[1:3,7:10]),ts);
    [imu, avp0] = avp2imu_ecef(avpe);
    avp = avptrans(avpe, 'e2g');  % err=avpcmpplot_polar(avp, inspure_grid(imu, avp0), 'g');
    if isfig
        imuplot(imu);
        insplot_polar(avp, 'g');
    end

