function [gpsPos, dxyz] = gpsnwpu(t0, T)
% 1Hz GNSS position obtained from real SPP sample @ NWPU within 48h.
%
% Prototype: [gpsPos, dxyz] = gpsnwpu(t0, T)
% Inputs: t0 - begin time in second, t0>1
%         T - time length in second, t0+T<48*3600
% Outputs: gpsPos - [lat, lon, height, t]
%          dxyz - xyz increment w.r.t. pos0 in meter
%          
% Example:
%   gps = gpsnwpu(13*3600, 3600*2);   gpsplot(gps);
%
% See also  gpssimu, dxyz2pos.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 24/12/2024
    load nwpu_gnss_1Hz48hur.mat;
    dxyz = double(dxyz)*20/2^7;
    gps = appendt(dxyz2pos(dxyz, pos0),1);
    gpsPos = gps(t0:t0+T,:);