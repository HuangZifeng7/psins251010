function dbu = vu2dbu(vu, isfig)
% Calculate up acc bias from up velocity.
%
% Prototype: dbu = vu2dbu(vu)
% Input: vu - up velocity, [vu,t]
% Outputs: dbu - up acc bias
%          isfig - figure flag
%
% See also  pp2vn, avplever.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 21/04/2025
global glv
    if nargin<2, isfig=1; end
    n = size(vu,2);
    if n>9, vu=vu(:,[6,end]);           % [avp,t]
    elseif n>6, vu=vu(:,[3,end]);  end  % [v/vp,t]
    dvu = diff(vu([1,end],:));
    dbu = dvu(1)/dvu(2);
    if isfig
        myfig,
        ax=plotyy(vu(:,2),vu(:,1),vu(:,2),vu(:,1)-(vu(:,2)-vu(1,2))*dbu);
        xyygo(ax,'t','VU',[labeldef('VU'),'_{res}']);
        title(sprintf('\\nabla_U=%.2f(ug)',dbu/glv.ug));
    end
