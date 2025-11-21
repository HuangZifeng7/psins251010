function out = bytecast(bytes,typ)
% Convert bytes to specific data type.
%
% Prototype: out = bytecast(bytes,typ)
% Input: bytes - uint8 type column bytes input
%        typ - int, uint, or float data output type.
% Outputs: out - specific output data array
%
% See also  byte2num, typecast.

% Copyright(c) 2009-2029, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 27/09/2025 
    bytes = uint8(bytes);
    [m,n] = size(bytes);
    if nargin<2
        if     n==2, typ='int16';
        elseif n==4, typ='int32';
        elseif n==8, typ='double'; end
    end
    bytes = reshape(bytes',1,m*n);
    out = typecast(bytes,typ)';

    % int32: bytes=bytes*[1;256;256^2;256^3]; idx=bytes>2^31 ?
