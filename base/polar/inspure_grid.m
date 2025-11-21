function avp = inspure_grid(imu, avp0, href, isfig)
% Process grid SINS pure inertial navigation with SIMU log data and
% using initial condition avp0 = [att0,vn0,pos0].
%
% Prototype: avp = inspure_grid(imu, avp0, href, isfig)
% Inputs: imu - SIMU data array
%         avp0 - initial parameters, avp0 = [att0,vn0,pos0]
%         href -    'v' - velocity fix-damping, =vn0
%                   'V' - vertical velocity fix-damping, =vn0(3)
%                   'p' - position fix-damping, =pos0
%                   'P' - position fix-damping, =pos0 & vertical velocity fix-damping, =vn0(3)
%                   'H' - height fix-damping, =pos0(3)
%                   'f' - height free.
%         isfig - figure on/off flag
% Output: avp - navigation results, avp = [attG,vG,pG,t]
%
% Example:
%   [imu, avp0, avp] = imupolar_grid(posset(89.9,0.10,100), 100, .1, 200);
%   avpg = inspure_grid(imu, avp0, 'V');
%   avpcmpplot_polar(avp, avpg, 'g');
%
% See also  insinit_grid, insupdate_grid, inspure, inspure_ecef, avpcmpplot_polar.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
global glv
    if nargin<3, href='f'; end
    vp_fix = href(1);
    [nn, ts, nts] = nnts(2, imu(:,end));
    ins = insinit_grid(avp0, ts);  vn0 = avp0(4:6); pos0 = avp0(7:9);  vG0 = ins.vG; pe0 = ins.pe;
    len = length(imu);    avp = zeros(fix(len/nn), 10);
    ki = timebar(nn, len, 'Grid-frame pure inertial navigation processing.');
    for k=1:nn:len-nn+1
        k1 = k+nn-1;
        wvm = imu(k:k1, 1:6);  t = imu(k1,end);
        ins = insupdate_grid(ins, wvm);
        if vp_fix=='v',      ins.vG = vG0;
        elseif vp_fix=='V',  ins.vG(3) = vG0(3);
        elseif vp_fix=='p',  ins.pe = pe0;
        elseif vp_fix=='P',  ins.pe = pe0;  ins.vG(3) = 0;
        elseif vp_fix=='H',  dh=ins.llh(3)-pos0(3);  ins.pe = ins.pe-dh*ins.uU;
        elseif vp_fix=='f',  ins.vG(3) = ins.vG(3);  % free, no need
        end
        avp(ki,:) = [ins.avp; t]';   %  avp(ki,1:3) = sum(imu(k:k1, 4:6))+ins.eth.gn'*0.01;
        ki = timebar;
    end
    if nargin<4, isfig=1; end
    if isfig==1,
        insplot_polar(avp,'g');
    end
