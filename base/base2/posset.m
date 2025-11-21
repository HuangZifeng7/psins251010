function pos = posset(pos0, lon, hgt, isdeg)
% Geographic position = [latitude; logititude; height] setting.
%
% Prototype: pos = posset(pos0, lon, hgt, isdeg)
% Input: pos0=[lat; lon; height], where lat and lon are always in arcdeg,
%             & height is in m.
%           or pos0=[lat; lon; hgt].
% Output: pos=[pos0(1)*arcdeg; pos0*arcdeg; hgt]
% 
% See also  avpset, llh, poserrset, initp.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 09/03/2014
global glv
    if nargin<4, isdeg=0; end 
    if     nargin>=3,  pos0 = [pos0; lon; hgt];
    elseif nargin==2,  pos0 = [pos0; lon; 0];  end  % pos = posset(lat, lon)
    if length(pos0)==1, pos0 = [pos0;0;0]; end      % pos = posset(lat)
    if     abs(pos0(1))>15959.999, isdeg=3;         % dms.sss
    elseif abs(pos0(1))>  159.999, isdeg=2;         % dm.mmm
    elseif abs(pos0(1))>     pi/2, isdeg=1; end     % d.ddd
    if isdeg==0
        pos = pos0;
    elseif isdeg==1
        pos = [pos0(1:2)*glv.deg; pos0(3)];
    elseif isdeg==2
        pos = [dm2r(pos0(1)); dm2r(pos0(2)); pos0(3)];
    elseif isdeg==3
        pos = [dms2r(pos0(1)); dms2r(pos0(2)); pos0(3)];
    end