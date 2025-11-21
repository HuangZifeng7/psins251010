function q = qmulBatch(q1, q2)
% Batch processing for quaternion multiplication: q = q1*q2.
% 
% Prototype: q = qmul(q1, q2)
% Inputs: q1, q2 - input quaternion
% Output: q - output quaternion ,such that q = q1*q2
%
% See also  qmul.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 14/05/2025
    isq3 = 0;
    if size(q1,2)==3, q1=q32q4(q1); isq3=1; end
    if size(q2,2)==3, q2=q32q4(q2); end
    if size(q1,1)~=size(q2,1), q2=q2(1:size(q1,1),:); end
    q = [ q1(:,1) .* q2(:,1) - q1(:,2) .* q2(:,2) - q1(:,3) .* q2(:,3) - q1(:,4) .* q2(:,4), ...
          q1(:,1) .* q2(:,2) + q1(:,2) .* q2(:,1) + q1(:,3) .* q2(:,4) - q1(:,4) .* q2(:,3), ...
          q1(:,1) .* q2(:,3) + q1(:,3) .* q2(:,1) + q1(:,4) .* q2(:,2) - q1(:,2) .* q2(:,4), ...
          q1(:,1) .* q2(:,4) + q1(:,4) .* q2(:,1) + q1(:,2) .* q2(:,3) - q1(:,3) .* q2(:,2) ];
    if isq3==1, q=q42q3(q); end