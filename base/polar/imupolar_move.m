function  [imu, avp0, avp] = imupolar_move(ap, ts, L0, isfig)
% Move low latitude trajectory to polar rigion.
%
% Prototype: [imu, avp0, avp] = imupolar_move(ap, ts, L0, isfig)
% Inputs: av - attENU,pos(lat,lon,h) in low latitude area
%         ts - sampling interval
%         L0 - mean latitude in (x+,0), default (0,0)
%         isfig - figure flag
% Outputs: imu - gyro & acc incremental outputs
%          avp0 - initial avp in ENU-frame
%          avp - avp series in grid-frame
%
% Example
%   load ap_flight.mat; ap=appendt([double(ap_a),ap_ll,double(ap_h)],ap_ts);  % insplot(ap,'ap');
%   [imu, avp0, avp] = imupolar_move(ap, 0.01, 89.9999*glv.deg);  % imuplot(imu);
%   avpg = inspure_grid(imu, avp0);  avpcmpplot_polar(avp, avpg);
%   avpn = inspure_enu(imu, avp0);   avpcmpplot_polar(avp, avptrans(avpn,'n2g'));
%
% See also  imupolar_rnd, imupolar_grid, yxyrnd, ap2avp, avptrans, avp2imu_ecef, aprot.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 14/12/2024
global glv
    [m,n] = size(ap);
    if nargin<4, isfig=1; end
    if nargin<3, L0=90*glv.deg; end
    if nargin<2, ts=0.01; end
    xyz = pos2dxyz(ap(:,end-3:end-1));  % north point to y+
    for k=1:3, xyz(:,k) = xyz(:,k)-mean(xyz(:,k)); end
    xyz(:,1) = xyz(:,1)+(90*glv.deg-L0)*glv.Rp;
    avp = [ap(:,1:2),ap(:,3)-pi/2,zeros(m,3),xyz(:,1:2),repmat(glv.Rp,m,1),ap(:,end)];  % -pi/2: grid-north point to x-
    llh=xyz2llhBatch(avp(:,7:9));  llh(:,3)=ap(:,end-1);  avp(:,7:9)=blh2xyzBatch(llh);  % insplot_polar(avp, 'g');
    avpe = avptrans(avp, 'g2e');  avpe = ap2avp_ecef(avpe(:,[1:3,7:10]),ts);
    [imu, avp0] = avp2imu_ecef(avpe);
    avp = avptrans(avpe, 'e2g');  % err=avpcmpplot_polar(avp, inspure_grid(imu, avp0), 'g');
    if isfig
        imuplot(imu);
        insplot_polar(avp, 'g');
    end

