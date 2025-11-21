function [dd, flen] = bffr(fname, mn, be)
% Binary file frame read, with fixed column length.
%
% Prototype: [dd, flen] = bffr(fname, mn, be)
% Inputs: mn - [m-row, n-column] to read
%         be - begin-end flag, 0 for begin, 1 for end, 2 for both begin&end
%                              >2 for every 'be' row read 1 row
% Outputs: dd - data array
%          flen - file length in byte
%
% See also  bffc, bffe, frm2num.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/09/2025 
    if length(mn)==1, mn=[1,mn]; end
    if nargin<3, be=0; end
    if be==-1, mn=[inf,mn(2)]; be=0; end  % read all
    fid = fopen(fname);
    dd = [];
    if be<=2
        if be==0 || be==2  % read first m-rows
            dd = fread(fid, [mn(2),mn(1)], 'uint8')';
        end
        if be==1 || be==2  % read last m-rows
            fseek(fid, -mn(1)*mn(2), 'eof');
            dd1 = fread(fid, [mn(2),mn(1)], 'uint8')';
            dd = [dd; dd1];
        end
        if nargout>1, fseek(fid, 0, 'eof'); flen = ftell(fid);  end
    else
        fseek(fid, 0, 'eof'); flen = ftell(fid); fseek(fid, 0, 'bof');
        rows = fix(flen/mn(2)/be);
        dd = zeros(rows,mn(2));
        for k=1:rows
            dd(k,:) = fread(fid, [mn(2),1], 'uint8')';
            fseek(fid, mn(2)*(be-1), 'cof');  % skip (be-1) rows
        end
    end
    fclose(fid);
