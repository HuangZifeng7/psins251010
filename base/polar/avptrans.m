function avp = avptrans(avp, str)
% AVP translation among ENU-frame, ECEF-frame and grid-frame.
%
% Prototype: avp = avptrans(avp, str)
% Inputs: avp - [att,vn,blh] in ENU-frame, [attECEF,vECEF,pe] in ECEF-frame
%               or [attG,vG,pe] in grid-frame
%         str - translation direction string, 'n2e','n2g','e2n','e2g','g2n','g2e'
%               'n' short for ENU-frame, 'g' for grid-frame, 'e' for ECEF-frame 
%               't' for transverse ENU-frame
% Output: avp - AVP output
%
% Example
%   avp0 = avpset([1;2;3], [4;5;6], [91;10;20]);
%   avp1 = avptrans(avp0,'n2n');
%   avp2 = [avp1(1:3)/glv.deg; avp1(4:6); avp1(7:8)/glv.deg; avp1(9)];
%
%   at = [0.1;0.2;1; 1;2;3; 0.1;0.2;1];
%   at = [0.1;0.2;0; 1;2;3; 0.1;1;1];
%   ae1 = avptrans(avptrans(at, 't2n'), 'n2e');
%   ae2 = avptrans(at, 't2e');  [ae1-ae2]
%
% See also  avp2ecefBatch, avp2gridBatch, blh2xyz, xyz2llh, pos2cng.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
global glv
    pos=avp(7:9); batch=0;
    if size(avp,2)>=9, pos=avp(1,7:9)'; batch=1; end
    if nargin<2, % default 'e2n' or 'n2e'
        if norm(pos)>glv.Re/2, str='e2n';
        else,                  str='n2e';  end
    end
    if batch==1
        switch str
            case {'n2e','e2n'}, avp=avp2ecefBatch(avp);
            case {'g2n','n2g'}, avp=avp2gridBatch(avp);
            case 'g2e', avp=avp2ecefBatch(avp2gridBatch(avp));  % g->n->e
            case 'e2g', avp=avp2gridBatch(avp2ecefBatch(avp));  % e->n->g
            case 'n2n', avp=avp2gridBatch(avp2gridBatch(avp));  % n->g->n
            case 'nen', avp=avp2ecefBatch(avp2ecefBatch(avp));  % n->e->n
            case 't2e',
                [Cen, avp(:,7:9)] = tpos2cen(avp(:,7:9));
                avp(:,4:6) = m3xv3(Cen,avp(:,4:6));
                avp(:,1:3) = m2attBatch(m3xm3(Cen, a2matBatch(avp(:,1:3))));
            case 'n2nc',  % n->g->n compare test
                avp1 = avp;
                avp = avptrans(avp, 'n2n');
                idx = avp1(:,7)>pi/2;
                avp1(idx,7)=pi-avp1(idx,7);  avp1(idx,[3,8])=avp1(idx,[3,8])+pi;  avp1(:,4:5)=-avp1(:,4:5); % NOTE!
                idx = avp1(:,3)>pi; avp1(idx,3)=avp1(idx,3)-2*pi;  idx = avp1(:,3)<-pi; avp1(idx,3)=avp1(idx,3)+2*pi;
                idx = avp1(:,8)>pi; avp1(idx,8)=avp1(idx,8)-2*pi;  idx = avp1(:,8)<-pi; avp1(idx,8)=avp1(idx,8)+2*pi;
%                 avp1 = avptrans(avp1, 'n2n');
                myfig
                subplot(321), plot([avp(:,1:2), avp1(:,1:2)]/glv.deg); xygo('pr');
                subplot(322), plot([avp(:,3), avp1(:,3)]/glv.deg); xygo('y');
                subplot(323), plot([avp(:,4:5), avp1(:,4:5)]); xygo('V');
                subplot(324), plot([avp(:,7), avp1(:,7)]/glv.deg); xygo('lat');
                subplot(325), plot([avp(:,8), avp1(:,8)]/glv.deg); xygo('lon');
        end
        return;
    end
    switch str
        case 'n2e',  %  avp-blh to ECEF-xyz
            xyz = blh2xyz(avp(7:9));
            Cen = pos2cen(avp(7:9));
            if length(avp)<9, avp(1:6) = [m2att(Cen*a2mat(avp(1:3))); xyz];
            else,             avp(1:9) = [m2att(Cen*a2mat(avp(1:3))); Cen*avp(4:6); xyz];  end            
        case 'n2g',  %  avp-blh to grid-xyz
            xyz = blh2xyz(avp(7:9));
            CGn = pos2cng(avp(7:9))';
            if length(avp)<9, avp(1:6) = [m2att(CGn*a2mat(avp(1:3))); xyz];
            else,             avp(1:9) = [m2att(CGn*a2mat(avp(1:3))); CGn*avp(4:6); xyz];  end
        case 'e2n',  %  ECEF-xyz to ENU-blh
            blh = xyz2llh(avp(7:9));
            Cne = pos2cen(blh)';
            if length(avp)<9, avp(1:6) = [m2att(Cne*a2mat(avp(1:3))); blh];
            else,             avp(1:9) = [m2att(Cne*a2mat(avp(1:3))); Cne*avp(4:6); blh];  end              
        case 't2n',  %  transverse-ENU to ENU
            [avp(7:9), Cnt] = tpos2pos(avp(7:9));
            avp(4:6) = Cnt*avp(4:6);
            avp(1:3) = m2att(Cnt*a2mat(avp(1:3)));
        case 'n2t',  %  ENU to transverse-ENU
            [avp(7:9), Ctn] = pos2tpos(avp(7:9));
            avp(4:6) = Ctn*avp(4:6);
            avp(1:3) = m2att(Ctn*a2mat(avp(1:3)));
        case 't2e',
             [Cen, avp(7:9)] = tpos2cen(avp(7:9));
             avp(4:6) = Cen*avp(4:6);
             avp(1:3) = m2att(Cen*a2mat(avp(1:3)));
        case 'e2g',
            avp = avptrans(avp, 'e2n');
            avp = avptrans(avp, 'n2g');
        case 'g2n',  % grid-xyz to ENU-blh
            blh = xyz2llh(avp(7:9));
            CnG = pos2cng(blh);
            if length(avp)<9, avp(1:6) = [m2att(CnG*a2mat(avp(1:3))); blh];
            else,             avp(1:9) = [m2att(CnG*a2mat(avp(1:3))); CnG*avp(4:6); blh];  end
        case 'g2e',
            avp = avptrans(avp, 'g2n');
            avp = avptrans(avp, 'n2e');
        case 'n2n',
            avp = avptrans(avp, 'n2g');
            avp = avptrans(avp, 'g2n');
    end
