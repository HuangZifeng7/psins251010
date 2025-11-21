function Cnb_FUR = CnbFUR(Cnb_RFU)
% Convert position matrix from ENU-frame to NUE-frame.
%
% Prototype: Cnb_FUR = CnbFUR(Cnb_RFU)
% Input: Cen_ENU - position matrix in ENU-frame
% Output: Cen_NUE - position matrix in NUE-frame
%
% See also  pos2cen, blh2xyz, xyz2blh.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 18/05/2025
    Cnb_FUR=[ Cnb_RFU(2,2),Cnb_RFU(2,3),Cnb_RFU(2,1);
              Cnb_RFU(3,2),Cnb_RFU(3,3),Cnb_RFU(3,1);
              Cnb_RFU(1,2),Cnb_RFU(1,3),Cnb_RFU(1,1) ];
