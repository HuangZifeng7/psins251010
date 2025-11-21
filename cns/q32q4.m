function q = q32q4(q)
% Translate q = [q1 q2 q3] to q = [q0 q1 q2 q3], for |q|=1 and q0>=0.
%
% Prototype: q = q32q4(q)
% Input: q - = [q1 q2 q3]
% Output: q - = [q0 q1 q2 q3], default for q0>=0
%
% Example
%    q=q32q4(randn(10,3)/3); [q, normv(q)]
%
% See also  q42q3

% Copyright(c) 2009-2021, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 09/11/2021
    if size(q,2)==1
        qq = q'*q;  if qq>1, q=q/sqrt(qq); qq=1.0; end
        qq1 = 1-qq;
        if qq1>=0,  q = [sqrt(qq1); q];
        else,       q = [0; q];      end
    else  % batch process
        sq = normv(q(:,1:3));
        idx = sq>=1;
        if sum(idx)>0, q(idx,1:3)=[q(idx,1)./sq(idx),q(idx,2)./sq(idx),q(idx,3)./sq(idx)]; sq(idx)=1.0; end
        q = [real(sqrt(1-sq.^2)), q];
    end