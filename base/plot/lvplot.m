function lvplot(lvx, lvy, lvz)
% Inner lever arm visualization plot 
%
% Example:
%   lvplot([-49.2; 3.5; 41.75], [-36.5;-21.2;41.75], [-34.5;5.5;68.45]);
%
% See also  imulever, avplever.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 31/03/2025
    if nargin<3, lvz = [0;0;0];  end
    if numel(lvx)==9, lvz=lvx(:,3); lvy=lvx(:,2); lvx=lvx(:,1);  end
    myfig
    plot3([0 lvx(1)], [0 lvx(2)], [0 lvx(3)], '-v', 'linewidth',2);  hold on;  grid on;  axis equal
    plot3([0 lvy(1)], [0 lvy(2)], [0 lvy(3)], '-^', 'linewidth',2);
    plot3([0 lvz(1)], [0 lvz(2)], [0 lvz(3)], '-d', 'linewidth',2);
    % ma = 1.1*max([0;lvx;lvy;lvz]);  mi = 1.1*min([0;lvx;lvy;lvz]);
    % xlim([-mi;ma]), ylim([-mi;ma]), zlim([-mi;ma]);
    xlabel('X'); ylabel('Y'); zlabel('Z'); 
    legend(['Ax',sprintf('%+9.3f ',lvx)], ['Ay',sprintf('%+9.3f ',lvy)], ['Az',sprintf('%+9.3f ',lvz)]);