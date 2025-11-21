function gylim(y0, y1, fign)
% In current SIMU figure (or fign), set all the gyro Y-axes limits to [y0, y1];
%
% Prototype: gylim(y0, y1, fign)
% Inputs: y0, y1 - ylim value.
%         fign - figure NO.
%
% See also  ylimall, xlimall, imuplot.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 01/10/2025
    if nargin<3, h=gcf; else, h=figure(fign); end
    if nargin<2,
        if nargin<1, y0=20; end  % 20dph
        if length(y0)<2, y1=y0;      y0=-y1;
        else             y1=y0(2);   y0=y0(1);   end
    end    
    for k=1:2:6
        subplot(3,2,k); ylim([y0,y1]);
    end
