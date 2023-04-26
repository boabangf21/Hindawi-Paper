function rayl = useraylmap(pos, raylmap, mapvec)
% DESCRIPTION rayl = useraylmap(pos, raylmap, mapvec)
% INPUT
%  pos --    Any complex matrix containing positions which may be interpreted as
%            absolute or relative positions in 2 dimensions.
%            any position is valid since wrapping into the map is done.
%            the standard deviation of the complex output rayl is 1 and the mean is 0.
%  mapvec -- two complex elements spanning the area of the map
% OUTPUT
%  rayl --   complex matrix of the same size as pos
%            the absolute value are rayleigh distributed 
% TRY       
%  [raylmap, mapvec] = creraylmap([1 i]*200);
%  plot(lin2db(abs(useraylmap(linspace(0,irandn,1000),raylmap,mapvec))))
%  mean(abs(useraylmap((0:0.01:100)*irand,raylmap,mapvec).^2))
% SEE ALSO  
%  creraylmap

% by Magnus Almgren 000229

% the same border must be used as when the map was created
border = 1;
siz = size(raylmap,1)-2*border; % size of side of map

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
rayl = interp2(mapi, mapi, raylmap.', pos2, pos1, 'cubic');

