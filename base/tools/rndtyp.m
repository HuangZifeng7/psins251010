function oldtyp = rndtyp(newtyp)
% Set random number generate type.
%
% Prototype: oldtyp = rndtyp(newtyp)
% Input: newtyp - new random type
% Output: oldtyp - old random type
%
% See also  rndgen, attrnd, posrnd.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 10/08/2025
global glv_rndtyp
    if ~exist('glv_rndtyp','var'), glv_rndtyp=2; end
    oldtyp = glv_rndtyp;
    glv_rndtyp = newtyp;