function ins = insinit_grid(avp0, ts)
% Grid SINS structure array initialization.
%
% Prototype: ins = insinit_grid(avp0, ts)
% Inputs: avp0 - initial avp0 = [att0; vn0; pos0] in n-frame
%         ts - SIMU sampling interval
% Output: ins - SINS structure array
%
% See also  insinit_ecef, insupdate_grid, insupdate, avpset.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/20224
global glv
    avp0 = avp0(:);
	ins = [];
	ins.ts = ts; ins.nts = 2*ts;
    [ins.pe, Cen] = blh2xyz(avp0(7:9));
    ins.CeG = pos2ceg(avp0(7:9));  CGn = ins.CeG'*Cen;
    ins.vG = CGn*avp0(4:6);
    ins.qGb = m2qua(CGn*a2mat(avp0(1:3)));
    ins.aG = zeros(3,1);

