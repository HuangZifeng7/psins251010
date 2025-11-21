function [pos, afa] = cew2pos(Cew)
% Convert geographic pos = [lat; lon; *] & wander azimuth to wander frame
%   transformation matrix Cew ( from Earth-frame to wander-frame ).
%
% Prototype: Cew = pos2cew(pos, afa)
% Inputs: pos - geographic position
%         afa - wander azimuth
% Output: Cew - transformation matrix from Earth-frame to wander-frame
%
% See also  pos2cew, pos2ceg, cew2pos.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 25/03/2025
    pos = [asin(Cew(3,3)); atan2(Cew(2,3),Cew(1,3)); 0];
    afa = atan2(-Cew(3,1),Cew(3,2));
    return;

    lat = 69*glv.deg;  afa = 0;
    Cew = pos2cew([lat;0;0]);
    wlon = 1*glv.dps;
    res = zeros(360,4);
    for k=1:360
        Cew = mupdt(Cew, rxyz(afa,'z')*[0;cos(lat)*wlon;0]);
        [pos, afa] = cew2pos(Cew);
        res(k,:) = [pos; afa]';
    end
    plotn(res/glv.deg);
