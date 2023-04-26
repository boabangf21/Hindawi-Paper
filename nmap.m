function nfold = nmap(rombvec, xyb, R_corrdist)
% DESCRIPTION nfold = nmap(rombvec, xyb, R_corrdist)
%  Decides what mapvector to use in crelogm.
% INPUT
%  xyb --         base stations positions
%  rombvec --     rombic vector
%  R_corrdist --  number of lognmap folding in rombvec
% OUTPUT
%  nfold --       The number of times that the map
%                 should fit into rombvec.
% SEE ALSO 
%  crelognmap, uselognmap, gcd1, gcd

% by Magnus Almgren 970310

% Find the minimum length vectors from site to site. 
xybb = xyb.' - xyb(round(size(xyb,2)/2));  
ind1 = find(xybb==0);
if length(ind1) > 0
 %Set all zero elements to max.
 xybb(ind1) = max(abs(xybb))*ones(size(ind1))+1;
end
ind2 = minind(abs(xybb),[],1); % This is now an index to a non zero difference.
if length(ind2) & max(xyb)>min(xyb)
 v = xybb(ind2)*[1; exp(i*pi/3)];  % A base pointing from site to site.
else
 nfold = 1; % Since there are only one site we do not have to fold the map.
 return
end
  
% ntimes is the number of v into rombvec
ntimes = [real(rombvec(1)).' imag(rombvec(1)).'] * inv([real(v) imag(v)]);
ntimes = round(ntimes);

% Split into rational numbers for nfold proposals between 1 and 100.
[a, b] = rat((1./(1:100))'*ntimes);  
ind3 = [1; find(max(b',[],1)'~=1)]; % the set of possible nfolds

% Select the one possible that is closest to the ratio between rombvec
% and number of hill * corrdist. 
nhill = 100;
ratio = abs(lin2db(min(abs(rombvec))./(ind3*sqrt(nhill)*R_corrdist)));

if all(ntimes)
 ntexp = repmat(abs(ntimes)', [1 length(ratio)]); % no 0
else
 ntexp = abs(ntimes(ntimes~=0)) * ones(1,length(ratio));
end

% Greatest divisors between ntimes and still possible nfolds.
divisor = gcd1([1:length(ratio); ntexp]);
possible = (divisor==1);
% We exclude the impossible cases by replacing values by Inf.
ratio(~possible) = ones(sum(~possible),1) * Inf;

% Minimum ratio of possible nfolds.
nfold = minind(ratio,[],1);
