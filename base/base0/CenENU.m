function Cen_ENU = CenENU(Cen_NUE)
% Convert position matrix from NUE-frame to ENU-frame.
%
% Prototype: Cen_ENU = CenENU(Cen_NUE)
% Input: Cen_ENU - position matrix in NUE-frame
% Output: Cen_NUE - position matrix in ENU-frame
%
% See also  CenNUE, pos2cen, blh2xyz, xyz2blh.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 24/06/2025
    Cen_ENU=Cen_NUE(:,[3,1,2]);
