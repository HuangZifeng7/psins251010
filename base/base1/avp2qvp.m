function  qvp = avp2qvp(avp, typ)
% Cubic spline interpolation of ap to generate avp with sampling time ts,
% where the velocity is the differentiation of position.

% Prototype: qvp = avp2qvp(avp, typ)
% Inputs: avp = [att, vn, pos, t] array
%         typ - =1 for launch vehicle, =0 for flight
% Output: qvp = [q1,q2,q3, vn, pos, t] array
%
% See also  qvp2avp, ap2avp.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 23/05/2025
    if nargin<2, typ=0; end
    qvp = avp;
    if typ==1
        if size(avp,2)==1
            qvp(1:3) = q42q3(a2qua1(avp(1:3)))';
        else
            qvp(:,1:3) = q42q3(a2qua1Batch(avp(:,1:3)))';
        end
    else
        if size(avp,2)==1
            qvp(1:3) = q42q3(a2qua(avp(1:3)))';
        else
            qvp(:,1:3) = q42q3(a2quaBatch(avp(:,1:3)))';
        end
    end
    