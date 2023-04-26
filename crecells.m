function [xyb, fib, rombvec] = crecells(cellradius, sps, km, lm, ncluster)
% DESCRIPTION [xyb, fib, rombvec] = crecells(cellradius, sps, km, lm, ncluster)
%  Creates basestations. Their positions in complex coordinates 
%  is stored in xyb. 
% INPUT
% cellradius --  Complex radius from the base to the center of the hexagonal cell.
% sps --         Sectors per site usually 1 (omni cells) or 
%                3 (120 degrees sector cells).
% km --          See below.
% lm --          Integer values giving the site-cluster-size according to
%                km^2+lm^2+km*lm. 
% nculster --    Square root of number of clusters created.
% OUTPUT
% xyb --         Base station positions in complex coordinates.
% fib --         Complex vectors from the bases to the center of the cell.
%                In the omni case fib vector elements are zero.
% rombvec --     Two vectors that span the area the system is created on.
% TRY            
%  plot(crecells(1000,3,2,3,2),'*') 
%  3 cells/site, 19 sites/cluster, 4 clusters 

% by Magnus Almgren 000505

[xys, clustervecs]= crecluster(km,lm);
xysites = flatten(xys,2); % 2nd dimension
ncl = 1:ncluster;
% Concatenate clusters in two dimensions.
cluv = flatten(mplus(ncl*clustervecs(1),ncl'*clustervecs(2)),4);
spsvec = flatten(0:sps-1,3); % 3rd dimension
fibd = cellradius*exp(i*2/3*pi*spsvec);  
[xyb, fib] = adjsiza(mplus(xysites,cluv), fibd);
xyb = cellradius*3*flatten(xyb,2); % base station positions 
xyb = xyb-mean(xyb); % adjust cells around origo
fib = flatten(fib,2);  % xyb & fib in second dimension
rombvec = clustervecs*ncluster*cellradius*3;
if sps == 1 % Modify for the omin case.
 fib = 0*fib;
 xyb = xyb/sqrt(3); % Adjust to correct cellradius.
 rombvec = rombvec/sqrt(3);
end 









