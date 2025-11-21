function r = rndgen(m,n)
% My special random number generate.
%
% Prototype: r = rndgen(m,n)
% Inputs: m,n - m-row, n-column
% Output: r - random number
%
% See also  rndtyp, attrnd, posrnd.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 10/08/2025
global glv_rndtyp
    if nargin<2, n=1; end
    if nargin<1, m=1; end
    switch glv_rndtyp
        case 0, r = zeros(m,n);
        case 1, r = ones(m,n);
        case 2, r = randn(m,n);
        case 3, r = rand(m,n);
        case 4, r = rand(m,n)*2-1;
        case 5, r = randi([0,1],m,n);
        case 6, r = randi([-1,1],m,n);
    end