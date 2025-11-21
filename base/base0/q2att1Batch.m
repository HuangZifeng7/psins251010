function att = q2att1Batch(qnb)
% Convert nX4 attitude quaternion to nX3 Euler attitude angles.
%
% Prototype: att = q2att1Batch(qnb)
% Input: qnb - attitude quaternion, nX4 array
% Output: att - Euler angles att=[pitch, roll, yaw] nX3 array in radians
%
% See also  q2att1, a2qua1, q2attBatch.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 14/05/2025

%     att = m2att1(q2mat(qnb));
    if size(qnb,2)==3, qnb=q32q4(qnb); end
    q11 = qnb(:,1).*qnb(:,1); q12 = qnb(:,1).*qnb(:,2); q13 = qnb(:,1).*qnb(:,3); q14 = qnb(:,1).*qnb(:,4); 
    q22 = qnb(:,2).*qnb(:,2); q23 = qnb(:,2).*qnb(:,3); q24 = qnb(:,2).*qnb(:,4);     
    q33 = qnb(:,3).*qnb(:,3); q34 = qnb(:,3).*qnb(:,4);  
    q44 = qnb(:,4).*qnb(:,4);
    C11=q11+q22-q33-q44;
    C21=2*(q23+q14);
    C31=2*(q24-q13); C32=2*(q34+q12); C33=q11-q22-q33+q44;
    att = [ atan2(C21,C11), atan2(C32,C33), -asin(C31) ];
    if size(qnb,2)>4, att(:,4) = qnb(:,end); end
    