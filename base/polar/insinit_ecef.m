function ins = insinit_ecef(avp0, ts)
% ECEF based SINS structure array initialization.
%
% Prototype: ins = insinit_ecef(avp0, ts)
% Inputs: avp0 - initial avp0 = [att0; vn0; pos0] in n-frame
%         ts - SIMU sampling interval
% Output: ins - SINS structure array
%
% See also  insupdate_ecef, insinit_grid, insinit_enu, insupdate, avpset.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2014
global glv
    avp0 = avp0(:);
	ins = [];
	ins.ts = ts; ins.nts = 2*ts;
    [ins.pe, Cen] = blh2xyz(avp0(7:9));
    ins.ve = Cen*avp0(4:6);
    ins.qeb = m2qua(Cen*a2mat(avp0(1:3)));
    ins.weie = [0;0;glv.wie];
    ins.ae = zeros(3,1);
    ins.sl = sin(avp0(7));
    ins.RNh = glv.Re/sqrt(1-glv.e2*ins.sl^2)+avp0(9);
    ins.ne = [ins.pe(1:2)/ins.RNh; ins.sl];

