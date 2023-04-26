function [lognmap, mapvec] = crelognmap(xyb, rombvec, corrdist)
% DESCRIPTION [lognmap, mapvec] = crelognmap(xyb, rombvec, corrdist); 
%  Generates a lognormal map used to generate lognormal fading
%  during a simulation. The upper border of the map has the same values 
%  as the lower border of the map. It is the same for the left and right borders
%  Within rombvec a certain number of maps are fitted to avoid resonace with base
%  station placement. The code is closely related to uselognmap. 
% INPUT
%  xyb --      Used to determine the size of the map.
%              If empty rombvec => mapvec. 
%  rombvec --  Used to determine the size of the map.
%  corrdist -- The distance until the correlation in the map has
%              decreased to 1/e.
% OUTPUT
%  lognmap --  A square matrix of real values. The mean of all values are 0
%              and the standard deviation is 1.
%  mapvec --   Two comlex vector elements that are used to move
%              the map to one mapsize.
% TRY 
%  mesh(crelognmap([],[1 i],0.1))
% SEE ALSO 
%  uselognmap, nmap

% by Magnus Almgren 000515

border = 1; % repeated border from the other side of map
siz = 2^8;  % size of side of map

mapvec = rombvec;
if ~isempty(xyb)
 mapvec = rombvec ./ nmap(rombvec, xyb, corrdist);
end

v1=linspace(-mapvec(1)/2, mapvec(1)/2, siz+1);
v2=linspace(-mapvec(2)/2, mapvec(2)/2, siz+1);

r = abs(mplus(fftshift(v2(1:siz)),fftshift(v1(1:siz)).'));
akf = exp(-r/corrdist);

% Calculate the map from white noise.
z = real(ifft2(fft2(randn(size(akf))).*sqrt(real(fft2(akf))))); 

% Extend it in each direction.  
lognmap = z([end-border+1:end 1:end 1:border],...
             [end-border+1:end 1:end 1:border]);



