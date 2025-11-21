function att = m2attBatch(Cnb)
% Convert direction cosine matrix(DCM) to Euler angles in batch processing.
%
% Prototype: att = m2attBatch(Cnb)
% Input: Cnb - =[C11, C12, C13, C21, C22, C23, C31, C32, C33, t]
% Output: att - =[pitch, roll, yaw, t]
%
% See also  a2mat, m2att, q2attBatch.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/11/2024
    att = q2attBatch(m2quaBatch(Cnb));
    if size(Cnb,2)==10, att(:,4)=Cnb(:,end); end
    