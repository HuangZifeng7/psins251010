function [xyz, Cen] = blh2xyzBatch(blh)
% Convert ECEF geographic coordinate [lat;lon;height] to
% Cartesian coordinate [x;y;z].  Note: BLH for Breite/Laenge/Hoehe?
%
% Prototype: [xyz, Cen] = blh2xyzBatch(blh)
% Input: blh - geographic coordinate blh=[lat;lon;height],
%               where lat & lon in radians and hegtht in meter
% Outputs: xyz - ECEF Cartesian coordinate vector, in meters
%          Cen - transformation matrix from Earth-frame to nav-frame
%
% See also  blh2xyz, xyz2blh, Dblh2Dxyz, pos2cen.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 29/11/2024
global glv
    B = blh(:,1); L = blh(:,2); H = blh(:,3);
    sB = sin(B); cB = cos(B); sL = sin(L); cL = cos(L);
    N = glv.Re./sqrt(1-glv.e2*sB.^2);
    X = (N+H).*cB.*cL;    Y = (N+H).*cB.*sL;    Z = (N*(1-glv.e2)+H).*sB;
    xyz = [X, Y, Z];
    if nargout==2
        Cen = [ -sL,  -sB.*cL,  cB.*cL, ...
                 cL,  -sB.*sL,  cB.*sL, ...
                 0*cB, cB,      sB ];
    end