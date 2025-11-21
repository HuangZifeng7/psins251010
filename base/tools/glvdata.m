function glvdata(row, clm)
% Save data to glv structure.
%
% Prototype: glvdata(row, clm)
% Inputs: row, clm - see the code
%
% Example
%   glvdata(10, 3);
%   plotn(glv.data);
%
% See also  glvf, psinslog.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 26/12/2024
global glv
    if nargin==1 % save data with only row
        if ~isfield(glv,'gdata'), glv.gdata=[]; glv.gcount=0; end  % init
        if size(glv.gdata,1)==glv.gcount, glv.gdata=[glv.gdata;zeros(size(glv.gdata))]; end  % double length
        glv.gcount = glv.gcount+1;
        glv.gdata(glv.gcount,:) = row(:)';
    elseif nargin==0  % delete residual
        glv.gdata(glv.gcount+1:end,:)=[];
    elseif nargin==2  % init row & column
        if size(row,1)>1, [row,clm]=size(row); end
        glv.gdata=zeros(row,clm); glv.gcount=0;
    end
    