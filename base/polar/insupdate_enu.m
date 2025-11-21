function ins = insupdate_enu(ins, imu)
% ENU-frame based SINS Updating Alogrithm.
%
% Prototype: ins = insupdate_enu(ins, imu)
% Inputs: ins - SINS structure array created by function 'insinit_enu'
%         imu - gyro & acc incremental sample(s)
% Output: ins - SINS structure array after updating
%
% See also  insinit_enu, insinit_ecef, insinit_grid, insupdate, insupdate_grid.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 11/12/2024
global glv
    nn = size(imu,1);
    nts = nn*ins.ts;  nts2 = nts/2;  ins.nts = nts;
    [phim, dvbm] = cnscl(imu,0);    % coning & sculling compensation
    ins.wib = phim/nts; ins.fb = dvbm/nts;
    %% earth & angular rate updating 
    vn01 = ins.vn+ins.an*nts2; pos01 = ins.pos+[vn01(2)/ins.RNh;vn01(1)/ins.clRNh;vn01(3)]*nts2;
    sl = sin(pos01(1));  cl = cos(pos01(1));
    tl = sl/cl;   sl2 = sl*sl;
    sq = 1-glv.e2*sl2;  sq2 = sqrt(sq);
    RMh = glv.Re*(1-glv.e2)/sq/sq2+pos01(3);
    ins.RNh = glv.Re/sq2+pos01(3);  ins.clRNh = cl*ins.RNh;
    ins.wnie = [0; glv.wie*cl; glv.wie*sl];
    ins.wnen = [-vn01(2)/RMh; vn01(1)/ins.RNh; vn01(1)/ins.RNh*tl];  % (*)
    ins.wnin = ins.wnie + ins.wnen;
    g = glv.g0*(1+5.2790414e-3*sl2+2.32718e-5*sl^4)-3.086e-6*pos01(3);
    gn01 = [0; 0; -g];
    %% (1)velocity updating
    gcc = -cross(2*ins.wnie+ins.wnen,vn01) + gn01;  % (**)
    ins.fn = qmulv(ins.qnb, ins.fb);
    ins.an = rotv(-ins.wnin*nts2, ins.fn) + gcc;
    vn1 = ins.vn + ins.an*nts;
    %% (2)position updating
    r = ins.clRNh;  r0=1000.00; % polar distance setting
    if r>r0  % non-polar region
        vn = (ins.vn+vn1)/2;  ins.vn = vn1;
        ins.pos = ins.pos + [vn(2)/ins.RNh;vn(1)/ins.clRNh;vn(3)]*nts;  % (***)
    else % near pole
        sl = sin(ins.pos(1));  cl = cos(ins.pos(1));   r = cl*ins.RNh;  % re-calcu r!
        ins.an = ins.fn+gn01;
        [dlon, ps, sy] = dllpolar(r, ins.vn, ins.an, nts);
%         ins.wnen = [(ps-r)/RMh; dlon*cl; dlon*sl]/nts;  ins.wnin = ins.wnie+ins.wnen;
        ins.wnen = [-sy/RMh; dlon*cl; dlon*sl]/nts;  ins.wnin = ins.wnie+ins.wnen;
        ins.vn = rotv(-(2*ins.wnie+ins.wnen)*nts, ins.vn+ins.an*nts);
        ins.pos = [sign(sl)*acos(ps/ins.RNh); ins.pos(2:3)+[dlon;ins.vn(3)*nts]];
    end
    if ins.pos(1)>pi/2, ins.pos(1)=pi-ins.pos(1); elseif ins.pos(1)<-pi/2, ins.pos(1)=-pi-ins.pos(1); end  % [-pi/2, pi/2]
    if ins.pos(2)>pi, ins.pos(2)=ins.pos(2)-2*pi; elseif ins.pos(2)<=-pi, ins.pos(2)=ins.pos(2)+2*pi; end % (-pi, pi]
    %% (3)attitude updating
    ins.qnb = qupdt2(ins.qnb, phim, ins.wnin*nts);
    ins.avp = [q2att(ins.qnb); ins.vn; ins.pos];
    ins.wUdlon = [ins.wnen(3); ins.vn(1)/ins.clRNh];
    