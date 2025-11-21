function magfind(app)
% Find mag dec/inc from https://www.magnetic-declination.com/.
%
% Prototype: magfind(app)
% Input: app - Windows application program case, msedge/chrome/firefox etc.
% Output: N/A
%
% See also  bdpoint, pos2bd, poscu.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 12/12/2024
    if nargin<1, app=1; end
    if app==1
        system('start msedge.exe https://www.magnetic-declination.com/');
    elseif app==2
        system('start chrome.exe https://www.magnetic-declination.com/');
    elseif app==3
        system('start firefox.exe https://www.magnetic-declination.com/');
    end