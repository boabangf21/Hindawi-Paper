 function [bter, kter, obkn] = terminate(b, k, obk, term)
% DESCRIPTION [bter, kter, obkn] = terminate(b, k, obk, term)
% INPUT
%  b --      Base index for each link. Size is #links by one. 
%  k --      Channel index for each link.  Size is #links by one.
%            The value NaN indicates not in use.
%  obk --    Binary matrix. Zero indicate occupied channel.
%  term --   A binary vector of the same size as b and k indicating
%            which links that are going to be released.
% OUTPUT
%  bter --   Has the same size as respectively b.
%  kter --   Nan in all places where calls where terminated. Has the same size as k.
%  obkn --   Reflects the new allocation in the system. Has the same size as obk.

% by Magnus Almgren 000505

if isempty(b)
 bter = zeros(0,1);
 kter = zeros(0,1);
 obkn = obk;
 break 
end
% NaN indicates not used any more
bter = b;
bter(term) = nan;
kter = k;
kter(term) = nan;
% release channels 
ind = index(b(term),k(term),size(obk));
obkn = obk;
obkn(ind) = obkn(ind) + 1; % Mark as available.
