function  avp = lciap2avp(ap, ts)
% Cubic spline interpolation of LCI-ap to generate avp with sampling time ts,
% where the velocity is the differentiation of position.

% Prototype: avp = lciap2avp(ap, ts)
% Inputs: ap = [att, x/y/z, t]
%         ts - sampling time interval
% Output: avp = [att, vx/y/z, x/y/z, t]
%
% See also  lciavp2imu, ap2avp, lcipure.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 02/05/2025
    ts0 = ap(2,7)-ap(1,7);
    if nargin<2,  ts = ts0;  end
    if size(ap,2)>7,  ap=ap(:,[1:3,end-3:end]);  end
    if ap(1,end)>0,  ap=[ap(1,1:3)-(ap(2,1:3)-ap(1,1:3))*ap(1,end)/ts0,0,0,0,0; ap];  end  % extend to t=0
    t = (ap(1,7):ts:ap(end,7))';
    avp = repmat(t,1,10);
    for k=1:3  % attitude spline interpolation
        avp(:,k) = spline(ap(:,7), ap(:,k), t);
    end
    for k=1:3
        pp = spline(ap(:,7), ap(:,3+k));
        avp(:,6+k) = ppval(pp, t); % position interpolation
        dpp = pp;
        for kk=1:pp.pieces
            dpp.coefs(kk,:) = [0, [3,2,1].*pp.coefs(kk,1:3)]; % coef differentiation
        end
        avp(:,3+k) = ppval(dpp, t); % velocity interpolation
    end
