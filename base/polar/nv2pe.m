function [pe, RN, RM] = nv2pe(ne, h)
% Convert normal vector & hieght to [px; py; pz] in ECEF frame.
%
% Prototype: [pe, RN] = nv2cen(ne, h)
% Inputs: ne - normal vector in ECEF frame
%         h - heigth
% Outputs: pe - = [px; py; pz] in ECEF frame
%          RN,RM - eastern,northern principal radius of curvature
%
% See also  nv2cen, pos2cen, cen2pos, pos2ceg, blh2xyz, xyz2blh, a2mat, pp2vn.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 21/04/2025
global glv
    if nargin<2, h=0; end
    I_e2sL2 = 1-glv.e2*ne(3)^2;  % 1-e^2*sinL^2
    RN = glv.Re/sqrt(I_e2sL2);
    pe = [(RN+h)*ne(1:2); (RN*(1-glv.e2)+h)*ne(3)];
    if nargout>2,  RM = RN*(1-glv.e2)/I_e2sL2;  end

