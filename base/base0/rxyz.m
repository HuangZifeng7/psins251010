function [C01, q01] = rxyz(ang, xyz)
% Rotation by x,y or z axis with angle 'ang'.
%
% Prototype: C01 = rxyz(ang, xyz)
% Inputs: ang - angle in rad
%         xyz - rotation axis 1/x/X, 2/y/Y or 3/z/Z
% Output: C01 - corresponding DCM = C^0_1
% 
% See also  rv2m, rv2q, a2mat, rotv.

% Copyright(c) 2009-2021, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 19/11/2021
    if nargin<2, xyz=3; end
	s = sin(ang); c = cos(ang);
    s2 = sin(ang/2); c2 = cos(ang/2);
    if size(ang,1)>1  % batch process
        O = zeros(length(ang),1);  I = ones(length(ang),1);
        switch xyz
            case {1, 'x', 'X'},
                C01 = [I O O, O c -s, O s c];
                q01 = [c2, s2,O,O];
            case {2, 'y', 'Y'}
                C01 = [c O s, O I O, -s O c];
                q01 = [c2, O,s2,O];
            case {3, 'z', 'Z'}
                C01 = [c -s O, s c O, O O I];
                q01 = [c2, O,O,s2];
        end
        return;
    end
    switch xyz
        case {1, 'x', 'X'},
            C01 = [1 0 0; 0 c -s; 0 s c];
            q01 = [c2; s2;0;0];
        case {2, 'y', 'Y'}
            C01 = [c 0 s; 0 1 0; -s 0 c];
            q01 = [c2; 0;s2;0];
        case {3, 'z', 'Z'}
            C01 = [c -s 0; s c 0; 0 0 1];
            q01 = [c2; 0;0;s2];
    end
