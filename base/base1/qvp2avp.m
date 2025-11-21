function  avp = qvp2avp(qvp, typ)
% Trans QVP t0 AVP.

% Prototype: qvp = avp2qvp(avp, typ)
% Inputs: qvp = [q1,q2,q3, vn, pos, t] array
%         typ - =1 for launch vehicle, =0 for flight
% Output: avp = [att, vn, pos, t] array with vn&pos untouched
%
% See also  avp2qvp, ap2avp.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 23/05/2025
    if nargin<2, typ=0; end
    avp = qvp;
    if typ==1
        if size(qvp,2)==1
            avp(1:3) = q2att1(q32q4(qvp(1:3)))';
        else
            avp(:,1:3) = q2att1Batch(q32q4(qvp(:,1:3)))';
        end
    else
        if size(avp,2)==1
            avp(1:3) = q2att(q42q3(qvp(1:3)))';
        else
            avp(:,1:3) = q2att1(q42q3(qvp(:,1:3)))';
        end
    end
