function [T, Dblh] = Dxyz2Dblh(blh, Dxyz)
% Convert perturbation error in ECEF Cartesian coordinate to geographic coordinate.
%
% Prototype: [T, Dblh] = Dxyz2Dblh(blh, Dxyz)
% Inputs: blh - geographic position [lat;lon;hgt]
%          Dxyz - position perturbation in Earth-frame in meters
% Outputs: T - transformation matrix from  Earth-frame to nav-frame, T=M_pv'*C^n_e
%         Dblh - position perturbation [dlat;dlon;dhgt]
%                where dlat,dlon in radians, dhgt in meter
%
% Example:
%   [T, Dblh] = Dxyz2Dblh(glv.pos0, [10;20;30]);
%   [T1, Dxyz] = Dblh2Dxyz(glv.pos0, Dblh);
%   [T*T1, Dxyz]
%
% See also  Dblh2Dxyz, blh2xyz, pos2cen.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 24/06/2025
global glv
    B = blh(1); L = blh(2); H = blh(3);
    sB = sin(B); cB = cos(B); sL = sin(L); cL = cos(L);
    N = glv.Re/sqrt(1-glv.e2*sB^2);
    a = N+H;  b = a-glv.e2*N;
    ba = b*cB^2+a*sB^2; aab = a-a*cB^2+b*cB^2;  bab = b+a*sB^2-b*sB^2;  acB=a*cB;
    T = [  -cL*sB/aab,   -sB*sL/bab,     cB/ba
              -sL/acB,       cL/acB,         0
          b*cB*cL/aab,  b*cB*sL/aab,   a*sB/ba ];
    if nargin==2
        Dblh = T*Dxyz;
    end

% syms B L a b
% sB=sin(B); cB=cos(B); sL=sin(L); cL=cos(L);
%     T = [ -a*sB*cL,         -a*cB*sL,  cB*cL;
%           -a*sB*sL,          a*cB*cL,  cB*sL;
%            b*cB,  0,         sB ];
% simplify(inv(T))
