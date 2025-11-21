function att = accmag2att(acc, mag)
% Using acc & mag to determine attitude in static base.
%
% Prototype: att = accmag2att(acc, mag)
% Inputs: acc - acc vector in b-frame, or acc&mag
%         mag - magnetic vector in b-frame, or [t1,t2]
% Output: att - Euler angles
%
% See also  sv2att, dv2atti, mv2atti, fb2atti, rv2q.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 23/01/2025
    if size(acc,2)>6  % att = accmag2att(accmag, [t1,t2])
        ag = mean(datacut(acc,mag(1),mag(2)),1);
        att = accmag2att(ag(end-6:end-4)', ag(end-3:end-1)');
        return;
    end
    [~, att, Cnb] = sv2atti([0;0;1], acc, 0);
    magn = Cnb*mag;
    att(3) = atan2(magn(1), magn(2));
