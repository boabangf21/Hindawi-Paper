function [xym, xyv] = cremob(nmob,rombvec)
% DESCRIPTION [xym, xyv] = cremob(nmob, rombvec)
%  Creates a number of uniformly distributed mobiles.
% INPUT
%  nmob --     Number of mobiles.
%  rombvec --  Two complex vectors spanning the rombic area containing cells.
% OUTPUT
%  xym --      Mobile positions, given as a complex column vector.
%  xyv --      Mobile velocity, given as a complex colum vector.
% TRY
%  [xym, xyv] = cremob(100,[1 i]), plot(xym, '*'), title('Mobile positions')
%  figure, plot(xyv, 'v'), title('Mobile velocity')

xym = rand(nmob, 2) * rombvec.';
xyv = irandn(size(xym))/sqrt(2)*14;