function [Mpv, Rx1, Ry1, tau1, slat, clat, slon, clon, ss, sc, cs, cc] = tpos2Mpv(tpos)
% Convert transverse tpos=[T_lat;T_lon;hgt] to matrix Mpv.
%
% Prototype: [Mpv, Rx1, Ry1, tau1, ...] = tpos2Mpv(tpos)
% Input: tpos - =[lat;lon;hgt] in transverse ENU frame
% Outputs: Mpv - velocity induced position change matrix
%          Rx1, Ry1, tau1 - curvature
%          slat, clat, slon, clon, ss, sc, cs, cc - sin/cos...
%
% See also  tpos2pos, tpos2cen, bl2sc.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 05/05/2025
global glv
    [slat, clat, slon, clon, ss, sc, cs, cc] = bl2sc(tpos);   cc2=cc*cc; sq1cc2=sqrt(1-cc2);
    if sq1cc2<1e-10
        cafa = 1; safa = 0;
    else
        cafa = -sc/sq1cc2;  safa = slon/sq1cc2;
    end
    RN = glv.Re/sqrt(1-glv.e2*cc2);  RM = RN*(1-glv.e2)/(1-glv.e2*cc2);
    RNh = RN+tpos(3);  RMh = RM+tpos(3);
    Rx1 = safa^2/RMh+cafa^2/RNh;  Ry1 = cafa^2/RMh+safa^2/RNh;  tau1 = safa*cafa*(1/RMh-1/RNh);
    Mpv = [-tau1,Ry1,0; [Rx1,-tau1]/clat,0; 0,0,1];