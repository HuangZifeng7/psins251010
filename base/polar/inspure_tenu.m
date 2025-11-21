function avp = inspure_tenu(imu, avp0, href, isfig)
% Process Transverse ENU-frame based SINS pure inertial navigation with 
% SIMU log data and using initial condition avp0 = [att0,vn0,pos0].
%
% Prototype: avp = inspure_tenu(imu, avp0, href, isfig)
% Inputs: imu - SIMU data array
%         avp0 - initial parameters, avp0 = [att0,vn0,pos0]
%         href -    'v' - velocity fix-damping, =vn0
%                   'V' - vertical velocity fix-damping, =vn0(3)
%                   'p' - position fix-damping, =pos0
%                   'P' - position fix-damping, =pos0 & vertical velocity fix-damping, =vn0(3)
%                   'H' - height fix-damping, =pos0(3)
%                   'f' - height free.
%         isfig - figure on/off flag
% Output: avp - navigation results, avp = [att,vn,pos,t]
%
% Example:
%   [imu, avp0, avp] = imupolar_ecef(posset(89.99,0,100), 10, .01, 200);
%   avpt = inspure_tenu(imuadderr(imu,imuerrset(0,0,0,000)), avp0, 'f');
%   avpcmpplot_polar(avp, avptrans(avpt,'t2e'), 'e');
%   avpn1 = inspure_enu(imu, avp0, 'V');  avpcmpplot_polar(avpn, avpn1, 'n');
%
% See also  insinit_tenu, insupdate_tenu, inspure, inspure_ecef, inspure_grid.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 06/02/2025
global glv
    if nargin<3, href='f'; end
    vp_fix = href(1);
    [nn, ts, nts] = nnts(2, imu(:,end));
    ins = insinit_tenu(avp0, ts);  vn0 = avp0(4:6); pos0 = avp0(7:9);
    len = length(imu);    avp = zeros(fix(len/nn), 10);
    ki = timebar(nn, len, 'ENU-frame pure inertial navigation processing.');
    for k=1:nn:len-nn+1
        k1 = k+nn-1;
        wvm = imu(k:k1, 1:6);  t = imu(k1,end);
        ins = insupdate_tenu(ins, wvm);
        if vp_fix=='v',      ins.vn = vn0;
        elseif vp_fix=='V',  ins.vn(3) = vn0(3);
        elseif vp_fix=='p',  ins.pos = pos0;
        elseif vp_fix=='P',  ins.pos = pos0;  ins.vn =vn0;
        elseif vp_fix=='H',  ins.pos(3)=pos0(3);
        elseif vp_fix=='f',  ins.vn(3) = ins.vn(3);  % free, no need
        end
        avp(ki,:) = [ins.avp; t]';
        ki = timebar;
    end
    if nargin<4, isfig=1; end
    if isfig==1,
        insplot(avp,'avp');  title('Tranverse ENU');
    end
