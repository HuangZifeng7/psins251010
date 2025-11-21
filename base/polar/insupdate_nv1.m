function ins = insupdate_nv1(ins, imu)
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
% 21/04/2025
global glv
    nn = size(imu,1);
    nts = nn*ins.ts;  nts2 = nts/2;  ins.nts = nts;
    [phim, dvbm] = cnscl(imu,0);    % coning & sculling compensation
    ins.wib = phim/nts; ins.fb = dvbm/nts;
    %% earth & angular rate updating 
    ve01 = ins.ve+ins.ae*nts2;  % extrapolation at t1/2
    ne01 = ins.ne - askew(ins.ne)*ins.M*ve01*nts;  ne01=ne01/norm(ne01);  sl=ne01(3);
    h01 = ins.h + ne01'*ve01*nts2;
    g = glv.g0*(1+5.2790414e-3*sl^2+2.32718e-5*sl^4)-3.086e-6*h01;
    ge01 = -g*ne01;
    %% (1)velocity updating
    gcc = -2*cross(ins.weie,ve01) + ge01;
    ins.fe = qmulv(ins.qeb, ins.fb);
    ins.ae = rotv(-ins.weie*nts2, ins.fe) + gcc;
    ve1 = ins.ve + ins.ae*nts;
    %% (2)position updating
    ve01 = (ins.ve+ve1)/2;  ins.ve = ve1;
    e2=glv.e2; ep2=glv.ep2; nx=ins.ne(1); ny=ins.ne(2); nz=ins.ne(3);
    RNh=glv.Re/sqrt(1-e2*nz^2)+ins.h; kL=(1-e2*nz^2)/(1-e2);
    ins.M = [-nx*ny*nz*ep2, -nz*(1+ep2*ny^2), ny*kL
             nz*(1+ep2*nx^2), nx*ny*nz*ep2, -nx*kL
             -ny, nx, 0]/RNh;
    ins.ne = ins.ne - askew(ins.ne)*ins.M*ve01*nts;  ins.ne=ins.ne/norm(ins.ne);  
    %ins.ne = ins.ne - askew(ins.ne)^2*ve01*nts/RNh;  ins.ne=ins.ne/norm(ins.ne);  
    ins.h = ins.h + ins.ne'*ve01*nts;
    %% (3)attitude updating
    ins.qeb = qupdt2(ins.qeb, phim, ins.weie*nts);
    if ins.qeb(1)<0, ins.qeb=-ins.qeb; end
    ins.avp = [q2att(ins.qeb); ins.ve; nv2pe(ins.ne,ins.h)];
