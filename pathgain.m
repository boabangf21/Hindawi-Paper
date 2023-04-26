function g = pathgain(xym, xyb, fib, rombvec, gainconst, alpha, sigma, raa, lobevector, lognmap, mapvec)
% DESCRIPTION g = pathgain(xym, xyb, fib, rombvec, gainconst, alpha, ...
%                          sigma, raa, lobevector, lognmap, mapvec)
%  Caluculates the gain in dB between all mobiles and all bases. 
% INPUT
%  xym --        Mobile coordinates, given as a complex column vector.
%  xyb --        Base station coordinates, given as a complex row vector.
%  fib --        Base station antenna vector, given as a row vector. Points
%                in the antenna direction.
%  rombvec --    Two complex vectors spanning the rombic area containing cells.
%  gainconst --  Gain at 1 meter distance.
%  alpha --      Distance attenuation coefficient.
%  sigma --      Standard deviation for the lognormal fading in dB.
%  raa --        Down link correlation (typical 0.5). This parameter determines
%                the correlation for the lognormal fading on the links
%                between the base stations and one mobile.
%  lobevector -- Antenna gain for all directions, a column vector 
%	              with 360 rows, one for each degree. 
%  lognmap --    Lognormal map used to create lognormal fading.
%  mapvec --     Two complex vectors that will determine the size of lognmap.
% OUTPUT
%  g --          Gain for all links in dB, given as a matrix 
%                one column for each base and one row for each mobil 
% TRY 
%  [lognmap, mapvec] = crelognmap([],[1 i]*10000,100);
%  g = pathgain(0,0,1000,[1 i]*10000,21,3.5,6,0.5,0,lognmap,mapvec)
%   Only one mobile and one base station.
%  xym = [rand(2)*[1;i]; rand(2)*[1;i]];
%  xyb = [1 i]*rand(2);
%  g = pathgain(xym,xyb,1000,[1 i]*10000,21,3.5,6,0.5,0,lognmap,mapvec)
%   Four mobiles and two base stations.
% SEE ALSO 
%  crelognmap, uselognmap

% by Magnus ALmgren 000510

if isempty(xym)
 g = zeros(0,size(xyb,2)); % Set to an empty matrix.
 return
end

% Wrap around in 9 fictive positions. 
bv = -1:1;
r = 0;  
if all(isfinite(rombvec)) % no rombvec, no wrap around
 r = flatten(mplus(rombvec(1)*bv.', rombvec(2)*bv),3);
end
wad = mplus(xym, -xyb, r);
 
% distance gain
% Make sure that the used distance is not less than 10 units (uma 95 08 02).
gatt = gainconst - alpha * 10 * log10(max(10,abs(wad)));

% Angle attenuation between direction to mobiles and cell main direction added.
gant=0;
fib(fib==0)=1;       % Point omnidirectional antennas west
bangle = angle(mdiv(wad,fib));
gant = lobevector(wrapind(bangle*180/pi,360));

% Select the position with highest gain of all wrap positions
g = max(gant+gatt,[],3);
clear gatt; % These may be big. 
clear gant;

% Lognormal fading with base station correlation added.
if sigma ~=0
 oseed = setseed(1); % Create a random base offset.
 maxdist = max(flatten(abs(mplus(xyb,-xyb.')))); % max distance 
 xyboffs = irand(size(xyb))*maxdist*10; % scramble base pos
 setseed(oseed);
 g = g + sigma * (...
  sqrt(1-raa) * uselognmap(mplus(-xym,xyboffs),lognmap,mapvec) + ...
  sqrt(raa)   * uselognmap(repmat(xym,[1 size(g,2)]),lognmap,mapvec)...
  );
end

















