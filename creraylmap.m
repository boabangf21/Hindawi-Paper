function [raylmap, mapvec] = creraylmap(rombvec)
% DESCRIPTION [raylmap, mapvec] = creraylmap(rombvec)
%  Creates a rayleigh fading map with complex samples.
%  The map is generated by filtering white noise through a bessel function.
%  This is made with a fourirer method that will make the map periodic.
%  This property is used by useraylmap.m to cover the simulated area with
%  a rayleighpattern much like a wall paper. mapvec is an integer fraction of 
%  rombvec. rombvec are two complex vectors that span the simulated area
% INPUT
%  rombvec --  The size of the area that is going to be covered by the map.
%              Two complex vectors spanning the simulated area.
% OUTPUT 
%  raylmap --  A small map with rayleigh samples that is used by 
%              the function useraylmap to create the actual
%              multipath fading.
%  mapvec --   Two complex elements spanning the area of the map.
% TRY
%  mesh(2*lin2db(abs(creraylmap([1 i]*1000))))
% SEE ALSO
%  useraylmap

% by Magnus Almgren 000229

ndips = 100; % average number of dips in the map.
lambda = 2e8/3e9; % wavelength at 2GHz
border = 1; % repeated border from the other side of map
siz = 2^7; % size of side of map in number of samples
mapvec = rombvec/ceil(min(abs(rombvec/lambda*2))/sqrt(ndips));
v = mprod(mapvec(:),linspace(-1/2, 1/2, siz));
r = abs(mplus(v(1,:).',v(2,:)))/lambda*2*pi;
htau = r.^-0.25.*besselj(0.25,r);

% calculate the map from white noise with convolution 
z = ifft2(fft2(irandn(size(htau))).*fft2(htau)); 
raylmap = z([1:end 1:border],[1:end 1:border]) ./std(z(:));  % normalize the map 
