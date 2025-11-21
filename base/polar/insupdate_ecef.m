function ins = insupdate_ecef(ins, imu)
% ECEF based SINS Updating Alogrithm.
% Ref: DuH, ECEF integrated navigation based on inertial and satellite navigation for polar region, OOT, 2022.
%      WangX．Global inertial navigation and integrated navigation method based on Earth coordinate frame, MCT, 2022.
%
% Prototype: ins = insupdate_ecef(ins, imu)
% Inputs: ins - SINS structure array created by function 'insinit_ecef'
%         imu - gyro & acc incremental sample(s)
% Output: ins - SINS structure array after updating
%
% See also  insinit_ecef, insinit_grid, insupdate, insupdate_enu, insupdate_grid.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
global glv
    nn = size(imu,1);
    nts = nn*ins.ts;  nts2 = nts/2;  ins.nts = nts;
    [phim, dvbm] = cnscl(imu,0);    % coning & sculling compensation
    ins.wib = phim/nts; ins.fb = dvbm/nts;
    %% earth & angular rate updating 
    ve01 = ins.ve+ins.ae*nts2; pe01 = ins.pe+ve01*nts2;  % extrapolation at t1/2
    [llh, RN, sl] = xyz2llh(pe01);
    g = glv.g0*(1+5.2790414e-3*sl^2+2.32718e-5*sl^4)-3.086e-6*llh(3);
    ge01 = -g*[pe01(1:2)/(RN+llh(3)); sl];
    %% (1)velocity updating
    gcc = -2*cross(ins.weie,ve01) + ge01;
    ins.fe = qmulv(ins.qeb, ins.fb);
    ins.ae = rotv(-ins.weie*nts2, ins.fe) + gcc;
    ve1 = ins.ve + ins.ae*nts;
    %% (2)position updating
    ins.pe = ins.pe + (ins.ve+ve1)*nts2;  
    ins.ve = ve1;
    %% (3)attitude updating
    ins.qeb = qupdt2(ins.qeb, phim, ins.weie*nts);
    if ins.qeb(1)<0, ins.qeb=-ins.qeb; end
    [ins.llh, ins.RN, ins.sl] = xyz2llh(ins.pe);
    ins.ne = [ins.pe(1:2)/(ins.RN+ins.llh(3)); ins.sl];  % normal vector in ECEF frame
    ins.avp = [q2att(ins.qeb); ins.ve; ins.pe];
    
%     Cen = pos2cen(ins.llh);
%     Cnb = Cen'*q2mat(ins.qeb);
%     ins.avp = [m2att(Cnb); Cen'*ins.ve; ins.llh];

