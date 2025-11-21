function ins = insupdate_nv(ins, imu)
% Normal Vector based SINS Updating Alogrithm.
% Ref: ZhangG, N-vector inertial navigation mechanization algorithm for transpolar aircraft, JCIT, 2017
%
% Prototype: ins = insupdate_nv(ins, imu)
% Inputs: ins - SINS structure array created by function 'insinit_ecef'
%         imu - gyro & acc incremental sample(s)
% Output: ins - SINS structure array after updating
%
% See also  insinit_nv, insinit_ecef, insinit_grid, insupdate, insupdate_ecef.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 05/02/2025
global glv
    nn = size(imu,1);
    nts = nn*ins.ts;  nts2 = nts/2;  ins.nts = nts;
    [phim, dvbm] = cnscl(imu,0);    % coning & sculling compensation
    ins.wib = phim/nts; ins.fb = dvbm/nts;
    %% earth & angular rate updating 
    ve01 = ins.ve+ins.ae*nts2;  % extrapolation at t1/2
    VeE = ins.Cen(:,1)'*ve01*ins.Cen(:,1);  VeN = ins.Cen(:,2)'*ve01*ins.Cen(:,2);
    weew = cros(ins.ne, VeE/(ins.RN+ins.h)+VeN/(ins.RM+ins.h)); 
    ne01 = ins.ne + cros(weew,ins.ne)*nts2;  ne01=ne01/norm(ne01);  sl=ne01(3);
    h01 = ins.h + ne01'*ve01*nts2;
    [Cen01, ~, RN01, RM01] = nv2cen(ne01, h01);
    g = glv.g0*(1+5.2790414e-3*sl^2+2.32718e-5*sl^4)-3.086e-6*h01;
    ge01 = -g*ne01;
    %% (1)velocity updating
    gcc = -2*cross(ins.weie,ve01) + ge01;
    ins.fe = qmulv(ins.qeb, ins.fb);
    ins.ae = rotv(-ins.weie*nts2, ins.fe) + gcc;
    ve1 = ins.ve + ins.ae*nts;
    %% (2)position updating
    ve01 = (ins.ve+ve1)/2;  ins.ve = ve1;
    VeE = Cen01(:,1)'*ve01*Cen01(:,1);  VeN = Cen01(:,2)'*ve01*Cen01(:,2);
    weew = cros(ne01, VeE/(RN01+h01)+VeN/(RM01+h01)); % angular rate: w^e_ew
    % ween1 = ve01'*Cen01(:,1)/(RN01+h01)*Cen01(:,2) - ve01'*Cen01(:,2)/(RM01+h01)*Cen01(:,1);
    ins.ne = ins.ne + cros(weew,ne01)*nts;  ins.ne=ins.ne/norm(ins.ne);  
    ins.h = ins.h + ins.ne'*ve01*nts;
    %% (3)attitude updating
    ins.qeb = qupdt2(ins.qeb, phim, ins.weie*nts);
    if ins.qeb(1)<0, ins.qeb=-ins.qeb; end
    [ins.Cen, ins.pe, ins.RN, ins.RM] = nv2cen(ins.ne, ins.h);
    ins.avp = [q2att(ins.qeb); ins.ve; ins.pe];
