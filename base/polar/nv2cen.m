function [Cen, pe, RN, RM] = nv2cen(ne, h)
% Convert normal vector & hieght to transformation matrix Cen ( from
%  Earth-frame to nav-frame £©.
%
% Prototype: [Cen, pe] = nv2cen(ne, h)
% Inputs: ne - normal vector in ECEF frame
%         h - heigth
% Outputs: Cen - transformation matrix from Earth-frame to nav-frame
%          pe - = [px; py; pz] in ECEF frame
%          RN,RM - eastern,northern principal radius of curvature
%
% See also  nv2pe, pos2cen, cen2pos, pos2ceg, blh2xyz, xyz2blh, a2mat, pp2vn.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 05/05/2025
global glv
    ne = ne/norm(ne);
    E = cros([0;0;1], ne);  % east vector
    if E'*E<1e-40, E(1)=1; end
    E = E/norm(E);  % normalization
    Cen = [ E, cros(ne,E), ne ];  % [ E, N, U ]
    if nargout>1
        if nargin<2, h=0; end
        I_e2sL2 = 1-glv.e2*ne(3)^2;  % 1-e^2*sinL^2
        RN = glv.Re/sqrt(I_e2sL2);  RM = RN*(1-glv.e2)/I_e2sL2;
        pe = [(RN+h)*ne(1:2); (RN*(1-glv.e2)+h)*ne(3)];
    end
