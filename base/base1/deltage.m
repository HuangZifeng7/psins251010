function [dge, dGe] = deltage(re, dre)
% Calculate the delta_ge from position vector re & error vector dre.
%
% Prototype: dge = deltage(re, dre)
% Inputs: re - position vector
%         dre - position error vector
% Outputs: dge - gravity error in ECEF-frame
%          dGe - universal gravitation in ECEF-frame
%
% See also  earth, gravj4.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 02/07/2025
global glv
    if re'*re<glv.Re
        pos = re;  [re, ~, RN] = blh2xyz(pos);
    else
        [pos, ~, RN] = xyz2blh(re);
    end
    ne = vnormlz(re);  sl = sin(pos(1));
    g = glv.g0*(1+5.2790414e-3*sl^2+2.32718e-5*sl^4)-3.086e-6*pos(3);
    D = 1./([RN;RN;RN*(1-glv.e2)]+pos(3));
    dge = -(g*diag(D)+g*(glv.e2*sl^2-1)*D.*ne*ne'-glv.beta2*(glv.beta*sl^2+1)*ne*ne')*dre +...
           (glv.g0*glv.e2*diag(D)-glv.beta2*glv.beta*eye(3))*sl*ne*dre(3);
    if nargin>1
        dGe = dge - askew([0;0;glv.wie])^2*dre;
    end
