function pos = cen2pos(Cen, h)
% Convert transformation matrix Cen to geographic pos = [lat; lon; hgt].
%
% Prototype: pos = cen2pos(Cen, h)
% Input: Cen - transformation matrix from Earth-frame to nav-frame
% Output: pos - geographic position
%
% See also  pos2cen, blh2xyz, xyz2blh, a2mat, pp2vn.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 05/05/2010
    if size(Cen,2)>=9  % batch process
        if nargin<2, h=Cen(:,end)*0; end
        pos = [asin(Cen(:,9)), atan2(Cen(:,6), Cen(:,3)), h];
        return;
    end
    if nargin<2, h=0; end
    pos = [asin(Cen(3,3)); atan2(Cen(2,3), Cen(1,3)); h];
