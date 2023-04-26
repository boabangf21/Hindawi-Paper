function [cdl, idl, sirdl]=transmitdl(b, k, pdl, g, noise)
% DESCRIPTION [cdl, idl, sirdl]=transmitdl(b, k, pdl, g, noise);
%  all variables in dB
% INPUT
%  b --      base station index
%  k --      channel index
%  pul --    transmitted power down link 
%  g --      gainmatrix
%  noise --  thermal noise floor [dBm]
% OUTPUT
%  cdl --    received carrier power down link [dBm]
%  idl --    interference power down link  [dBm]
%  sirdl --  signal to interference ratio down link [dB]
% TRY 
%  [cdl, idl, sirdl]=transmitdl([1 2]', [1 1]', [0 0]', [-2 -10; -10 -3], -inf)

% by Magnus Almgren 000510

cdl = nans(size(k)); % prepare for result
idl = nans(size(k));
if all(isnan(k))
 sirdl = cdl - idl; 
 return
end
tl = isfinite(k); %these are true links
m = (1:length(k))';
ind = index(m(tl), b(tl),size(g,1)); % index into g possibly containing nans
gcl = g(ind); % vector, gain for connected links
cdl(tl) = pdl(tl)+g(ind); % carrier
sdbk = db2lin(g) * ... % a b by k matrix with received power
 sparse(b(tl),k(tl),db2lin(pdl(tl)),size(g,2),max(k)) + db2lin(noise); 
% interference per link
idl(tl) = lin2db(max(0,sdbk(index(m(tl), k(tl), size(g,1))) - db2lin(cdl(tl)))); %    Received interference, i = ss - c
sirdl = cdl - idl; % C/I
