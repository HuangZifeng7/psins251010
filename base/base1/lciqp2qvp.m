function  qvp = lciqp2qvp(qp, ts)
% Cubic spline interpolation of LCI-qp to generate qvp with sampling time ts,
% where the velocity is the differentiation of position.

% Prototype: qvp = lciap2avp(qp, ts)
% Inputs: qp = [q1/q2/q3, x/y/z, t]
%         ts - sampling time interval
% Output: qvp = [q1/2/3, vx/y/z, x/y/z, t]
%
% See also  lciqvp2imu, lciap2avp, ap2avp.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 12/05/2025
    ts0 = qp(2,7)-qp(1,7);
    if nargin<2,  ts = ts0;  end
    if size(qp,2)>7,  qp=qp(:,[1:3,end-3:end]);  end
    if qp(1,end)>0,  qp=[qp(1,1:3)-(qp(2,1:3)-qp(1,1:3))*qp(1,end)/ts0,0,0,0,0; qp];  end  % extend to t=0
    t = (qp(1,7):ts:qp(end,7))';
    qvp = repmat(t,1,11);
    quat = [sqrt(1-normv(qp(:,1:3))),qp(:,1:3)];  % asume q0>=0
    for k=1:4  % quat spline interpolation
        qvp(:,k) = spline(qp(:,7), quat(:,k), t);
    end
    qvp(:,2:4) = qvp(:,2:4)./normv(qvp(:,1:4));
    qvp(:,1) = [];  % delete q0
    for k=1:3
        pp = spline(qp(:,7), qp(:,3+k));
        qvp(:,6+k) = ppval(pp, t); % position interpolation
        dpp = pp;
        for kk=1:pp.pieces
            dpp.coefs(kk,:) = [0, [3,2,1].*pp.coefs(kk,1:3)]; % coef differentiation
        end
        qvp(:,3+k) = ppval(dpp, t); % velocity interpolation
    end
