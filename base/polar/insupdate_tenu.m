function ins = insupdate_tenu(ins, imu)
% Transverse geographic coordinate ENU-frame based SINS Updating Alogrithm.
% Ref: XuX, Inertial navigation algorithm in polar regions based on transerse geographic coordinate system, JHUST, 2014
%
% Prototype: ins = insupdate_tenu(ins, imu)
% Inputs: ins - SINS structure array created by function 'insinit_enu'
%         imu - gyro & acc incremental sample(s)
% Output: ins - SINS structure array after updating
%
% See also  insinit_tenu, insinit_enu, insinit_ecef, insinit_grid, insupdate, insupdate_grid.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 06/02/2025
global glv
    nn = size(imu,1);
    nts = nn*ins.ts;  nts2 = nts/2;  ins.nts = nts;
    [phim, dvbm] = cnscl(imu,0);    % coning & sculling compensation
    ins.wib = phim/nts; ins.fb = dvbm/nts;
    %% earth & angular rate updating 
    vn01 = ins.vn+ins.an*nts2; pos01 = ins.pos+ins.Mpv*vn01*nts2;
    [ins.Mpv, Rx1, Ry1, tau1, slat, clat, slon, clon, ss, sc, cs, cc] = tpos2Mpv(pos01);
    ins.wnie = glv.wie*[-slon; -sc; cc];
    ins.wnen = [tau1,-Ry1; Rx1,-tau1; slat/clat*[Rx1,-tau1]]*vn01(1:2);
    ins.wnin = ins.wnie + ins.wnen;
    g = glv.g0*(1+5.2790414e-3*cc^2+2.32718e-5*cc^4)-3.086e-6*pos01(3);
    gn01 = [0; 0; -g];
    %% (1)velocity updating
    gcc = -cross(2*ins.wnie+ins.wnen,vn01) + gn01;
    ins.fn = qmulv(ins.qnb, ins.fb);
    ins.an = rotv(-ins.wnin*nts2, ins.fn) + gcc;
    vn1 = ins.vn + ins.an*nts;
    %% (2)position updating
    ins.pos = ins.pos + ins.Mpv*(ins.vn+vn1)*nts2;
    ins.vn = vn1;
    %% (3)attitude updating
    ins.qnb = qupdt2(ins.qnb, phim, ins.wnin*nts);
    ins.avp = [q2att(ins.qnb); ins.vn; ins.pos];
    