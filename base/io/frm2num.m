function num = frm2num(frm, frmstr)
% Trans data frame to double number array
%
% Prototype: num = frm2num(frm, frmstr)
% Inputs: frm - uint8 data array
%         frmstr - frame string discription: i/u/f+8/16/32/64bits, k+1/2/3...bytes
% Outputs: num - double output data array
%
% Example
%   num = frm2num(randn(10,150), 'u8,i32,k5,f32,f64');
%
% See also  binfrmfile, byte2num, bffr, bffe.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 28/09/2025 
    frm = uint8(frm);
    for k=length(frmstr):-1:1
        if frmstr(k)==' '||frmstr(k)==','||frmstr(k)==';', frmstr(k)=[]; end
    end
    iuf='iufk'; idx=[];
    for k=1:length(iuf), idx=union(idx,strfind(frmstr,iuf(k))); end
    num = [];  clm = 1;
    for k=1:length(idx)
        typ = frmstr(idx(k));  bits = sscanf(frmstr(idx(k)+1:end),'%d');
        if typ=='i', typ='int';
        elseif typ=='u', typ='uint';
        elseif typ=='f'&&bits==32, typ='single';
        elseif typ=='f'&&bits==64, typ='double';  end
        if strcmp(typ,'int')||strcmp(typ,'uint'), typ=[typ,sprintf('%d',bits)]; end
        if typ(1)=='k'
            num = [num, zeros(size(frm,1),1)];  clm = clm+bits;  % 'k' for skip 'bits' bytes
        else
            num = [num, double(bytecast(frm(:,clm:(clm+bits/8-1)),typ))];  clm = clm + bits/8;
        end
    end
