function str = ptitle(str1, num1, str2, num2, str3, num3, str4, num4)
% PSINS title setting for plot
%
% Prototype: str = ptitle(str1, num1, str2, num2, str3, num3, str4, num4)
% Inputs: stri, num1 - specific string and numerical setting.
%
% See also  xygo, labeldef, xyygo, myfig, xlimall, nextlinestyle.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 23/08/2025
    str = [labeldef(str1),' = ',sprintf('\\rm%.3f ', num1)];
    if nargin>2
        str = [str, '; ', labeldef(str2),' = ',sprintf('\\rm%.3f ', num2)];
    end
    if nargin>4
        str = [str, '; ', labeldef(str3),' = ',sprintf('\\rm%.3f ', num3)];
    end
    if nargin>6
        str = [str, '; ', labeldef(str4),' = ',sprintf('\\rm%.3f ', num4)];
    end
    title(str);


