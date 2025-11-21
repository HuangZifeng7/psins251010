function vec = vnormlz(vec)
% Vector normalization, so ||vec||=1.
%
% Prototype: vec = vnormlz(vec)
% Input: vec - input vector whose norm may not be 1
% Output: vec - input vector whose norm equals 1
%
% See also  mnormlz, qnormlz, normv.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 23/09/2012
    [m, n] = size(vec);
    if m>1 && n>1  % 04/12/2024, for row-vectors
        nm = normv(vec);
        for k=1:n, vec(:,k) = vec(:,k)./nm; end
        return;
    end
    vec = vec/norm(vec);