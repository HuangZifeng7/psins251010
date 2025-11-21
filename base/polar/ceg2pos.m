function pos = ceg2pos(CeG, h)
% Convert transformation matrix CeG (from Earth-frame to grid-frame £©to
% geographic pos = [lat; lon; hgt].
%
% Prototype: pos = ceg2pos(CeG, h)
% Inputs: CeG - transformation matrix from Earth-frame to grid-frame
%         h - height
% Output: pos - geographic position [lat,lon,hgt]
%
% See also  pos2ceg, pos2cen, cen2pos, blh2xyz, xyz2blh.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
    if nargin<2, h=0; end
    pos = [asin(CeG(3,3)); atan2(CeG(2,3),CeG(1,3)); h];