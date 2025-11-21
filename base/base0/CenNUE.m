function Cen_NUE = CenNUE(Cen_ENU)
% Convert position matrix from ENU-frame to NUE-frame.
%
% Prototype: Cen_NUE = CenNUE(Cen_ENU)
% Input: Cen_ENU - position matrix in ENU-frame
% Output: Cen_NUE - position matrix in NUE-frame
%
% See also  CenENU, pos2cen, blh2xyz, xyz2blh.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 02/05/2025
    Cen_NUE=[ Cen_ENU(1,2),Cen_ENU(1,3),Cen_ENU(1,1);  % =Cen_ENU(:,[2,3,1]);
              Cen_ENU(2,2),Cen_ENU(2,3),Cen_ENU(2,1);
              Cen_ENU(3,2),Cen_ENU(3,3),Cen_ENU(3,1) ];
    % Cen = [ -slon,  -slat*clon,  clat*clon
    %          clon,  -slat*slon,  clat*slon
    %          0,      clat,       slat      ];
