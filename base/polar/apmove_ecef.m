function  ap = apmove_ecef(ap0, pos0)
% Move original trajectory ap0 to a specific place, whose first point
% is at pos0 and initial yaw is yaw0.
%
% Prototype: ap = apmove(ap0, pos0, yaw0)
% Inputs: ap0 - original att&pos parameters
%         pos0 - new initial position
%         yaw0 - new initial yaw
% Output: ap - att&pos after linear translation % angular rotation 
%
% See also  apscale, ap2imu, avp2imu, posstatic, trjsimu, insupdate.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 16/03/2014
    blh0 = xyz2blhBatch(ap0(:,4:6));
    xyz0 = blh2xyz(pos0);
    rv = angle3d(ap0(1,4:6)'/norm(ap0(1,4:6)), xyz0/norm(xyz0));
    Cee = rv2m(rv);
    Ceb0 = a2matBatch(ap0(:,1:3));
    Ceb = [ Cee(1,1)*Ceb0(:,1)+Cee(1,2)*Ceb0(:,4)+Cee(1,3)*Ceb0(:,7), ...
            Cee(1,1)*Ceb0(:,2)+Cee(1,2)*Ceb0(:,5)+Cee(1,3)*Ceb0(:,8), ...
            Cee(1,1)*Ceb0(:,3)+Cee(1,2)*Ceb0(:,6)+Cee(1,3)*Ceb0(:,9), ...
            Cee(2,1)*Ceb0(:,1)+Cee(2,2)*Ceb0(:,4)+Cee(2,3)*Ceb0(:,7), ...
            Cee(2,1)*Ceb0(:,2)+Cee(2,2)*Ceb0(:,5)+Cee(2,3)*Ceb0(:,8), ...
            Cee(2,1)*Ceb0(:,3)+Cee(2,2)*Ceb0(:,6)+Cee(2,3)*Ceb0(:,9), ...
            Cee(3,1)*Ceb0(:,1)+Cee(3,2)*Ceb0(:,4)+Cee(3,3)*Ceb0(:,7), ...
            Cee(3,1)*Ceb0(:,2)+Cee(3,2)*Ceb0(:,5)+Cee(3,3)*Ceb0(:,8), ...
            Cee(3,1)*Ceb0(:,3)+Cee(3,2)*Ceb0(:,6)+Cee(3,3)*Ceb0(:,9) ];
    blh = xyz2blhBatch(ap0(:,4:6)*Cee');
    ap = [m2attBatch(Ceb), blh2xyzBatch([blh(:,1:2),blh0(:,3)]), ap0(:,7)];
    
  