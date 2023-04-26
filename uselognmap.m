function logn = uselognmap(pos, lognmap, mapvec)
% DESCRIPTION logn = uselognmap(pos, lognmap, mapvec)
%  Generates correlated lognormal values that are position dependent.
%  This routine is closely related to crelognmap.m
% INPUT
%  pos --        any matrix with complex positions.
%                values are folded back into the map when outside
%  lognmap --    Global variable containing the lognormal map.
%  mapvec --     two complex vector elements describing the size of the map
% OUTPUT
%  logn --       Value of lognormal fading. Same size as the "pos"
%                input argument
% SEE ALSO 
%  crelognmap

% by Magnus Almgren 960508, revised 991101

% the same border must be used as when the map was create
border = 1; % extension for cubic interpolation
siz = size(lognmap,1)-2*border; % size of side of map

% the inverse base of v (v(1)*base => [1 0] and v(2)*base => [0 1] 
base = inv([real(mapvec.') imag(mapvec.')]);

% the portions of v(1) and v(2) for each position element
rpos = [real(pos(:)) imag(pos(:))];

% wrap it  back into the interval {0,1}
mpos = mod(rpos*base,1);
pos1 = reshape(mpos(:,1), size(pos));
pos2 = reshape(mpos(:,2), size(pos));

% vectors used by interp2 for table lookup 
mapi = (-1:siz)/siz;

% make the interpolation
logn = interp2(mapi, mapi, lognmap.', pos2, pos1, 'cubic');
