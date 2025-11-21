function ins = insinit_enu(avp0, ts)
% ENU-frame based SINS structure array initialization.
%
% Prototype: ins = insinit_enu(avp0, ts)
% Inputs: avp0 - initial avp0 = [att0; vn0; pos0] in n-frame
%         ts - SIMU sampling interval
% Output: ins - SINS structure array
%
% See also  insupdate_enu, insinit_ecef, insinit_grid, insupdate, avpset.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 11/12/2024
global glv
    avp0 = avp0(:);
	ins = [];
	ins.ts = ts;  ins.iter = 0;
    ins.qnb = a2qua(avp0(1:3)); ins.vn = avp0(4:6); ins.pos = avp0(7:9);
    ins.an = zeros(3,1);
    sl = sin(ins.pos(1));  cl = cos(ins.pos(1));
    if abs(cl)<1/glv.Rp, cl=1/glv.Rp; end
    sq = 1-glv.e2*sl^2;  sq2 = sqrt(sq);
    ins.RNh = glv.Re/sq2+ins.pos(3);  ins.clRNh = cl*ins.RNh;
