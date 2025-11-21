function phi = pn2phistd(pn, lti, alignT, intervalT)
% Calculating yaw misalign angle STD from pure SINS velocity error.
%
% Prototype: phi = vn2phistd(vn, lti, alignT, intervalT)
% Inputs: vn - pure SINS velocity error, in most case for static base
%         lti - latitude
%         alignT - align time lenght
%         intervalT - interval time
% Output: phi - misalignment angles
%
% See also  vn2phistd.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 17/01/2025
global glv
    if nargin<3, alignT=300; end
    if nargin<4, intervalT=alignT/5; end
    if intervalT<1, intervalT=1; elseif intervalT>alignT, intervalT=alignT; end
    if length(lti)>1, lti=lti(1); end  % lti = pos;
    t0 = pn(1,end); pn(:,end)=pn(:,end)-t0;
    n = fix((pn(end,end)-alignT)/intervalT)+2;
    phi = zeros(n,4);
    for k=1:n
        t1 = (k-1)*intervalT; t2 = t1+alignT;
        if t2>pn(end,end)+1, phi(k:end,:)=[]; break; end
        [~, ~, phik] = pn2phiu(datacut(pn,t1,t2), lti, 0);
        phi(k,:) = [phik; t1]';
    end
    pn(:,end) = pn(:,end)+t0;
    phi(:,end) = phi(:,end)+t0;
    myfig,
    subplot(221), plot(pn(:,end), pn(:,1:2)); xygo('dP'); xm=xlim;
    title(sprintf('alignT=%.1fs, intervalT=%.1fs',alignT,intervalT));
    subplot(223), plot(phi(:,end), phi(:,3)/glv.min, '-o'); xygo('phiU'); xlim(xm);
    title(sprintf('std=%.2f^\\prime',std(phi(:,3))/glv.min));
    subplot(222), plot(phi(:,end), phi(:,1)/glv.sec, '-o'); xygo('phiE'); xlim(xm);
    title(sprintf('std=%.2f^\\prime^\\prime',std(phi(:,1))/glv.sec));
    subplot(224), plot(phi(:,end), phi(:,2)/glv.sec, '-o'); xygo('phiN'); xlim(xm);
    title(sprintf('std=%.2f^\\prime^\\prime',std(phi(:,2))/glv.sec));

