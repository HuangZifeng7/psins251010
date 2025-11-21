function qrb = trjattrot(qrb, rotstr, w, ts, st, rep)
% Body frame attitude trajectory simulation.
%
% Prototype: qrb = trjattrot(qrb, rotstr, w, ts, st, rep)
% Inputs: qrb - initial quaternion
%         rotstr - rotation string
%         w - rotation angular rate
%         ts - sampling interval
%         st - static time before&after each rotation
%         rep - rotation repeat
% Output: qrb - trajectory quaternion array, [q0,q1,q2,q3, t]
%
% Example
%    qrb = trjattrot([], 'y10x20y-30', 10*glv.dps, 0.01, 10, 10);
%    qrb = trjattrot(qrb, 'z179.9z-179.9', 10*glv.dps, 0.01, 30, 10);
%    insplot(qrb, 'q');
%
% See also  rot2ang, attrottt, pfa2qua, trjSimu.

% Copyright(c) 2009-2025, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 13/05/2025
    if nargin<6, rep=1; end
    if nargin<5, st=10; end
    if length(st)==1; st=[st;st]; end
    if isempty(qrb), qrb=[1;0;0;0;0]'; end
    if size(qrb,2)==1, qrb=[qrb;0]'; end  % [q0;q1;q2;q3] , default t=0
    t = qrb(end,end);
    rotstr = upper(rotstr);
    x = strfind(rotstr,'X');  y = strfind(rotstr,'Y');  z = strfind(rotstr,'Z');
    xyz = union(union(x,y),z);
    for r=1:rep
        for k=1:length(xyz)
            ang = sscanf(rotstr(xyz(k)+1:end),'%f')*pi/180;
            if abs(ang)>1*pi/180/60
                ang = rot2ang(ang, w, ts, 10);
            end
            ang = [repmat(ang(1,1),fix(st(1)/ts),1); ang(:,1); repmat(ang(end,1),fix(st(2)/ts),1)];
            len = length(ang);
            q0 = qrb(end,:)';  qrbk = zeros(len,5);
            for m=1:len
                [C, q] = rxyz(ang(m),rotstr(xyz(k)));
                qk = qmul(q0, q);  t = t+ts;
                qrbk(m,:) = [qk; t]';
            end
            qrb = [qrb; qrbk];
        end
    end
    % qrb = [q42q3(qrb),(1:length(qrb))'*ts];
