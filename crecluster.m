function [xys, clv] = crecluster(km,lm)
% DESCRIPTION [xys, clv] = crecluster(km,lm)
%  Creates one cluster containing km^2 + lm^2 + km*lm sites.
%  These terms are handled as three parallellograms (two rhombs) 
%  and then combined to a cluster. There is no equal channels 
%  inside a cluster.
% INPUT
%  km --   See below.
%  lm --   Integer. At least one of km and lm is >0. 
% OUTPUT
%  xys --  The positions of the sites in complex coordinates.
%  clv --  Two vectors used to repeat clusters next 
%          to each other. The vectors tells the displcement in complex 
%          coordinates from the cluster to its nearest copies.
% TRY 
%  [xys, clustervecs] = crecluster(3,2), plot(xys,'*')
%  figure, plot(mplus(xys.',[0 -clustervecs clustervecs]),'*')

% by Magnus Almgren 990716

e = exp(i*2/3*pi*(0:2));   % vectors used to span the cluster
a = e*exp(i*pi/6)/sqrt(3); % vectors used to span the cluster

kv = 0:km-1; 
lv = 0:lm-1;
[m1 m2] = meshgrid(kv*e(1),kv*e(2));
skk = flatten(m1+m2)+a(1);
[m1 m2] = meshgrid(lv*e(1),lv*e(3));
sll = flatten(m1+m2)+a(3);
[m1 m2] = meshgrid(kv*e(2),lv*e(3));
slk = flatten(m1+m2)+a(2);
clv = (km*e(1)-lm*e(3))*[1 e(2)];
xys = [skk(:).' sll(:).' slk(:).'];





