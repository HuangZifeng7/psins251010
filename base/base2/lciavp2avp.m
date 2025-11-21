function avpn = lciavp2avp(avpi, pos0, A0, typ)
% Transform LCI AVP to ENU AVP.
%
% Prototype: avpn = lciavp2avp(avpi, pos0, A0)
% Inputs: avpi - LCI-frame AVP =[pitch,roll,yaw, vx,vy,vz, x,y,z, t],
%         pos0 - initial launch position [lat;lon;hgt]
%         A0 - firing/shooting angle, NOTE: counter-clockwise to be +
%         typ - avpi type define
%               typ==0:  RFU yaw-pitch-roll for att, 
%                        x-rightward, y-forward, z-upward for vel&pos
%               typ==1:  FUR pitch-yaw-roll for att, 
%                        x-forward, y-upward, z-rightward for vel&pos
% Output: avpn - ENU AVP out =[pitch,roll,yaw, ve,vn,vu, lat,lon,hgt, t]
%
% See also  lciavp2imu, atttrans, q2att1, axxx2a.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 01/05/2025
global glv
    if nargin<4, typ=0; end
    if typ==0
        Cn0a = a2mat([0;0;A0]);
        [re0,Ce0n0] = blh2xyz(pos0);  Ce0a = Ce0n0*Cn0a;  ra0 = Ce0a'*re0;
    else
        Cn0a = a2mat([0;A0;0]);
        [re0,Ce0n0] = blh2xyz(pos0);  Ce0n0=CenNUE(Ce0n0);  Ce0a = Ce0n0*Cn0a;  ra0 = Ce0a'*re0;
    end
    avpn = avpi;
    len = length(avpi);
    if len>10000, timebar(1, len, 'Trans LCI AVP to ENU AVP.'); end
    for k=1:len
        att = avpi(k,1:3)';  va = avpi(k,4:6)';  pa = avpi(k,7:9)';  t = avpi(k,10)';
        Ceka = rv2m([0;0;-glv.wie*t])*Ce0a;
        if typ==0
            re = Ceka*(ra0+pa);  [pos,Ceknk] = xyz2blh(re);  Cnka = Ceknk'*Ceka;
            avpn(k,:) = [m2att(Cnka*a2mat(att)); Ceknk'*(cross([0;0;-glv.wie],re)+Ceka*va); pos; t];
        else
            re = Ceka*(ra0+pa);  [pos,Ceknk] = xyz2blh(re);  Ceknk=CenNUE(Ceknk); Cnka = Ceknk'*Ceka;
            [~,Cab] = a2qua1(att);
            Cnkb = Cnka*Cab;
            Cnkb = [ Cnkb(3,3),Cnkb(3,1),Cnkb(3,2);  % C^n_b FUR->RFU
                     Cnkb(1,3),Cnkb(1,1),Cnkb(1,2);
                     Cnkb(2,3),Cnkb(2,1),Cnkb(2,2) ];
            vnk = Ceknk'*(cross([0;0;-glv.wie],re)+Ceka*va);
            avpn(k,:) = [m2att(Cnkb); vnk([3,1,2]); pos; t];
        end
        if len>10000, timebar; end
    end

