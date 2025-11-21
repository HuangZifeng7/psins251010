function ins = insinit_tenu(avp0, ts)
% Transverse ENU-frame based SINS structure array initialization.
%
% Prototype: ins = insinit_tenu(avp0, ts)
% Inputs: avp0 - initial avp0 = [att0; vn0; pos0] in n-frame
%         ts - SIMU sampling interval
% Output: ins - SINS structure array
%
% See also  insupdate_tenu, insinit_enu, insinit_ecef, insinit_grid, insupdate, avpset.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 06/02/2025
global glv
    avp0 = avp0(:);
	ins = [];
	ins.ts = ts;
    avp0 = avptrans(avp0, 'n2t');
    ins.qnb = a2qua(avp0(1:3)); ins.vn = avp0(4:6); ins.pos = avp0(7:9);
    ins.an = zeros(3,1);
    ins.Mpv = tpos2Mpv(ins.pos);
