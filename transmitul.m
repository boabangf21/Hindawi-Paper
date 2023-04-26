function [cul, iul, sirul]=transmitul(b, k, pul, g, noise)
% DESCRIPTION [cul, iul, sirul]=transmitul(b, k, pul, g, noise)
%  Calculates SIR uplink.
% INPUT
%  b --      base station index
%  k --      channel index
%  pul --    transmitted power up link 
%  g --      gainmatrix
%  noise --  thermal noise floor [dBm]
% OUTPUT
%  cul --    received carrier power up link [dBm]
%  iul --    interference power up link  [dBm]
%  sirul --  signal to interference ratio up link [dB]
% TRY 
%  [cdl, idl, sirdl]=transmitdl([1 2]', [1 1]', [0 0]', [-2 -10; -10 -3], -inf)

% by Magnus Almgren 000517

cul = nans(size(k)); % prepare for result
iul = nans(size(k));
if all(isnan(k))
 sirul = cul - iul; 
 return
end
tl = isfinite(k); %these are true links
m = (1:length(k))';
ind = index(m(tl), b(tl),size(g,1)); % index into g possibly containing NaNs
gcl = g(ind); % vector, gain for connected links
cul(tl) = pul(tl)+g(ind); % carrier
subk = db2lin(g') * ... % a b by k matrix with received power
 sparse(m(tl),k(tl),db2lin(pul(tl)),size(g,1),max(k)) + db2lin(noise); 
% interference per link
iul(tl) = lin2db(max(0,flatten(subk(index(b(tl), k(tl), size(g,2)))) - db2lin(cul(tl)))); % Received interference, i = ss - c
sirul = cul - iul; % C/I

