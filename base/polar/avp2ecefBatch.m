function avp = avp2ecefBatch(avp, n2e)
% Batch processing AVP translation, from ENU-frame to ECEF-frame,
% or reverse.
%
% Prototype: avp = avp2ecefBatch(avp, n2e)
% Input: avp - [att,vn,blh] in ENU-frame, if norm(blh)>Re/2 for in ECEF-frame
%        n2e - ==1 from ENU-frame to ECEF-frame, ==0 reverse
% Output: avp - [qeb(2:4),ve,pe] in ECEF-frame, or in ENU-frame
%
% See also  avp2grid, blh2xyz, xyz2llh.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
global glv
    if nargin<2, n2e=1; end
    if norm(avp(1,7:9))>glv.Re/2, n2e=0; end
    if n2e==0  % ECEF-xyz to avp-blh
        for k=1:size(avp,1)
            blh = xyz2llh(avp(k,7:9));
            Cne = pos2cen(blh)';
            avp(k,1:9) = [m2att(Cne*a2mat(avp(k,1:3))); Cne*avp(k,4:6)'; blh]';
        end
        return;
    end
%     idx = avp(:,7)>pi/2;
%     avp(idx,7)=pi-avp(idx,7);  avp(idx,[3,8])=avp(idx,[3,8])-pi;  avp(:,4:5)=-avp(:,4:5); % NOTE!
    B = avp(:,7); L = avp(:,8); H = avp(:,9);
    sB = sin(B); cB = cos(B); sL = sin(L); cL = cos(L);
    N = glv.Re./sqrt(1-glv.e2*sB.^2);
    X = (N+H).*cB.*cL;    Y = (N+H).*cB.*sL;    Z = (N*(1-glv.e2)+H).*sB;
    pe = [X, Y, Z];
    Cen = [ -sL,  -sB.*cL,  cB.*cL,...
             cL,  -sB.*sL,  cB.*sL,...
             0*B,  cB,      sB ];
    ve = [ Cen(:,1).*avp(:,4)+Cen(:,2).*avp(:,5)+Cen(:,3).*avp(:,6),...
           Cen(:,4).*avp(:,4)+Cen(:,5).*avp(:,5)+Cen(:,6).*avp(:,6),...
           Cen(:,7).*avp(:,4)+Cen(:,8).*avp(:,5)+Cen(:,9).*avp(:,6) ];
    Cnb = a2matBatch(avp(:,1:3));
    Ceb = [ Cen(:,1).*Cnb(:,1)+Cen(:,2).*Cnb(:,4)+Cen(:,3).*Cnb(:,7), ...
            Cen(:,1).*Cnb(:,2)+Cen(:,2).*Cnb(:,5)+Cen(:,3).*Cnb(:,8), ...
            Cen(:,1).*Cnb(:,3)+Cen(:,2).*Cnb(:,6)+Cen(:,3).*Cnb(:,9), ...
            Cen(:,4).*Cnb(:,1)+Cen(:,5).*Cnb(:,4)+Cen(:,6).*Cnb(:,7), ...
            Cen(:,4).*Cnb(:,2)+Cen(:,5).*Cnb(:,5)+Cen(:,6).*Cnb(:,8), ...
            Cen(:,4).*Cnb(:,3)+Cen(:,5).*Cnb(:,6)+Cen(:,6).*Cnb(:,9), ...
            Cen(:,7).*Cnb(:,1)+Cen(:,8).*Cnb(:,4)+Cen(:,9).*Cnb(:,7), ...
            Cen(:,7).*Cnb(:,2)+Cen(:,8).*Cnb(:,5)+Cen(:,9).*Cnb(:,8), ...
            Cen(:,7).*Cnb(:,3)+Cen(:,8).*Cnb(:,6)+Cen(:,9).*Cnb(:,9) ];
    avp(:,1:9) = [q2attBatch(m2quaBatch(Ceb)), ve, pe];
