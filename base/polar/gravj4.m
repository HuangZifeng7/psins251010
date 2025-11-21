function [G, ge, gn, err] = gravj4(re, Ji)
% Calculate gravity by using J0~J6 model.
%
% Prototype: [f, ge, gn, err] = gravj4(re, Ji)
% Inputs: re - [x; y; z] in ECEF frame
%         Ji - using Ji model flag, i=0,2,4,6, default=4
%              =0 for J0 model, =2 for J2, =3 for J4, =6 for J6
% Outputs: G - universal gravitation in i-frame
%          ge - gravity in ECEF frame
%          gn - gravity in ENU frame
%          err - error between Somigliana model & Ji model
% 
% See also  grav, egmwgs84, earth, xyz2blh.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 21/04/2025
global glv
    if nargin<2, Ji=4; end
    r = sqrt(re'*re);  z = re(3);
    gmr = glv.GM/r^3;
    f0 = -gmr*re;  f2=zeros(3,1); f4=f2; f6=f2;
    if Ji>=2
        zr = z/r;  zr2 = zr*zr;   Rr2 = glv.Re^2/r^2;
        f2 = 3/2*glv.J2*gmr*Rr2 * ((5*zr2-1)*re-[0;0;2*z]);
    end
    if Ji>=4
        zr4 = zr2*zr2;  Rr4 = Rr2*Rr2;
        f4 = 5/8*glv.J4*gmr*Rr4 * (3*(21*zr4-14*zr2+1)*re-[0;0;4*z*(7*zr2-3)]);
    end
    if Ji>=6
        zr6 = zr4*zr2;  Rr6 = Rr4*Rr2;
        f6 = 7/16*glv.J6*gmr*Rr6 * ((429*zr6-495*zr4+135*zr2-5)*re-[0;0;6*z*(33*zr4-30*zr2+5)]);
    end
    G = f0+f2+f4+f6;
    if nargout>1
        wr = glv.wie^2*[re(1:2);0];
        ge = G+wr;
    end
    if nargout>2
        [pos, Cen] = xyz2blh(re);
        gn = Cen'*ge;
        eth = earth(pos);
        err = gn-eth.gn;
        % C = zeros(11);  S = C;
        % C(1:2:11,1) = [ 1;                              -0.108262982131*10^-2/sqrt(5); % WGS-84 coefficients
        %                 0.237091120053*10^-5/sqrt(9);   -0.608346498882*10^-8/sqrt(13); 
        %                 0.142681087920*10^-10/sqrt(17); -0.121439275882*10^-13/sqrt(21) ];
        % gn1 = egm(glv.GM, glv.Re, glv.wie, glv.f, C, S, pos(1), pos(2), pos(3));
        % err = gn-gn1;
    end
    return;

   %% Example
   glvs; errk=[]; % glv.wie=0;
   for k=1:89.999, [G,ge,gn,err]=gravj4(blh2xyz([k*glv.deg;0*glv.deg;0]),6); errk(k,:)=[err;k]'; end
   myfig, plot(errk(:,end), errk(:,1:3)/glv.ug);  xygo('L / \circ', 'err / ug');

    %% diff verify - debug
    dxyz=3; dx=0; dy=0; dz=0; 
    if dxyz==1, dx=1e-1; elseif dxyz==2, dy=1e-1; else, dz=1e-1; end
    pos = [89*glv.deg;0;000];  xyz = blh2xyz(pos);
    x=xyz(1); y=xyz(2); z=xyz(3); r=norm([x;y;z]);  zr=z/r;
    V00 = glv.GM/r;
    V20 = -glv.J2*glv.GM*glv.Re^2/2*(3*z^2/r^5-1/r^3);
    V40 = -glv.J4*glv.GM*glv.Re^4/8*(35*z^4/r^9-30*z^2/r^7+3/r^5);
    V60 = -glv.J6*glv.GM*glv.Re^6/16*(231*z^6/r^13-315*z^4/r^11+105*z^2/r^9-5/r^7);
    x1=x+dx; y1=y+dy; z1=z+dz;  r1=norm([x1;y1;z1]);
    V01 = glv.GM/r1;
    V21 = -glv.J2*glv.GM*glv.Re^2/2*(3*z1^2/r1^5-1/r1^3);
    V41 = -glv.J4*glv.GM*glv.Re^4/8*(35*z1^4/r1^9-30*z1^2/r1^7+3/r1^5);
    V61 = -glv.J6*glv.GM*glv.Re^6/16*(231*z1^6/r1^13-315*z1^4/r1^11+105*z1^2/r1^9-5/r1^7);
    dV0 = (V01-V00)./[dx;dy;dz];
    dV2 = (V21-V20)./[dx;dy;dz];
    dV4 = (V41-V40)./[dx;dy;dz];
    dV6 = (V61-V60)./[dx;dy;dz];
    f0 = -glv.GM/r^3*[x;y;z];
    f2 = [  3*glv.J2*glv.GM*glv.Re^2/2/r^5*(5*zr^2-1)*x;
            3*glv.J2*glv.GM*glv.Re^2/2/r^5*(5*zr^2-1)*y;
            3*glv.J2*glv.GM*glv.Re^2/2/r^5*(5*zr^2-3)*z  ];
    f4 = [ 15*glv.J4*glv.GM*glv.Re^4/8/r^7*(21*zr^4-14*zr^2+1)*x;
           15*glv.J4*glv.GM*glv.Re^4/8/r^7*(21*zr^4-14*zr^2+1)*y;
           15*glv.J4*glv.GM*glv.Re^4/8/r^7*(21*zr^4-14*zr^2+1)*z - ...
           20*glv.J4*glv.GM*glv.Re^4/8/r^7*(7*zr^2-3)*z ];
    f6 = [ 7*glv.J6*glv.GM*glv.Re^6/16/r^9*(429*zr^6-495*zr^4+135*zr^2-5)*x;
           7*glv.J6*glv.GM*glv.Re^6/16/r^9*(429*zr^6-495*zr^4+135*zr^2-5)*y;
           7*glv.J6*glv.GM*glv.Re^6/16/r^9*(429*zr^6-495*zr^4+135*zr^2-5)*z - ...
          42*glv.J6*glv.GM*glv.Re^6/16/r^9*(33*zr^4-30*zr^2+5)*z ];
    f00 = [dV0, f0]/glv.ug
    f22 = [dV2, f2]/glv.ug
    f44 = [dV4, f4]/glv.ug
    f66 = [dV6, f6]/glv.ug
