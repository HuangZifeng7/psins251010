function [phiu0, vnfit, phi] = pn2phiu(pn, lti, isfig)
% Calculating yaw misalign angles from pure SINS position error(in meter).
%
% Prototype: [phiu0, vnfit, phi] = vn2phiu(vn, lti, isfig)
% Inputs: pn - pure SINS position error, pn(:,1)=east, pn(:,2)=north (in meter)
%         lti - latitude
%         isfig - figure flag
% Output: phiu0 - misalignment between calculating navigation frame and real
%               navigation frame
%         vnfit - polyfit for velocity
%         phi - [phiE0, phiN0, phiU0];
%
% See also  vn2phiu.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 20/01/2025
    if nargin<3, isfig=0; end
    t0 = 2*pn(1,end)-pn(2,end); t = pn(:,end)-t0;  vn(:,4)=t(1:10:end);
    p1 = polyfit(t,pn(:,1),3); vn(:,1) = polyval(polyder(p1),vn(:,4));
    p2 = polyfit(t,pn(:,2),3); vn(:,2) = polyval(polyder(p2),vn(:,4));
    p3 = polyfit(t,pn(:,3),3); vn(:,3) = polyval(polyder(p3),vn(:,4));
    vn(:,4) = vn(:,4)+t0;
    [phiu0, vnfit, phi] = vn2phiu(vn, lti, isfig);
    return;
    
     myfig, subplot(211), plot(vn(:,end), vn(:,1:3));  legend('E','N','U');
     subplot(212), plot(pn(:,end), pn(:,1:3)-pn(1,1:3), vn(:,end),cumsum(vn(:,1:3))*diff(vn(1:2,end)));
