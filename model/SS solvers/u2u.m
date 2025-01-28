function [vv1,vv2,vv3] = u2u(p,uu1,uu2,uu3)
%% function [vv1,vv2,vv3] = u2u(p,uu1,uu2,uu3)
% where: [uu] == [n;m;I];
% 
% Syntax
% 
% [vv1,vv2,vv3] = u2u(p,uu1) returns three P-by-N matrices n, m and I where
% P is the number of space steps, and N is the number of columns of uu1.
% Each column of uu1 corresponds to a different state. 
% 
% [vv1] = u2u(p,uu1,uu2,uu3) returns the state matrix uu1. Each column of
% uu1 is a state vector corresponding to n,m and I vertically concatenated 

switch nargin
    case 2
        vv1 = uu1(1:p.P,:);
        vv2 = uu1(p.P+1:2*p.P,:);
        vv3 = uu1(2*p.P+1:3*p.P,:);
    case 4
        vv1 = [uu1;uu2;uu3];
    otherwise
        error('wrong number of arguments');
end