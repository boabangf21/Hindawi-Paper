function [bass,kass,offern] = assign(b,k,g,offer,homargin)
% DESCRIPTION [bass,kass,offern] = assign(b,k,g,offer,homargin)
%  Assigns base and channel to new calls. If possible
%  first base is selected and then if possible a channel.
%  The assignment order is random, thus the order within b and k 
%  has no impact in a blocking situation.
% INPUT
% b --          Base index for each link. Size is #links by one.
%               Elements of b that is equal to NaN are going to be assigned
%               a base station index.
% k --          Channel index for each link. Size is #links by one.
%               The value NaN indicates no channel selected yet. 
% g --          Gain matrix with gain values in dB. Size is #links by #bases.
% offer --      Binary matrix. Zero indicate occupied channel. 
%               Size is #bases by #channels.
% homargin --   Standard deviation of the noise added to the
%               gain matrix when base selction is made.
% OUTPUT
% bass --       As b, but possibly with NaNs removed by an index
%               to indicate the assigned base station.             
% kass --       As k, but possibly with NaNs removed by an index
%               to indicate the assigned channel. If blockin occours,
%               i.e. no channel available on the selected base, then the
%               corresponding element will be kept as nan.
% offern --     Same size as offer, but with assigned base station 
%               channel combination marked as occupied (zero).
% TRY 
%  [b,k,offer]=assign(nan,nan,-100,1,0)
%  [b,k,offer]=assign(nans(3,1),nans(3,1),10*rand(3,2),ones(2,4),3)
%  3 new users in a sytem with 2 bases and with the same 4 channels in 
%  both bases. The gain matrix is just log-normal in order to indicate some 
%  gain variation. Reduce the number of channels to one and see what happens.

%  by Magnus Almgren 000505

if isempty(b) % the empty case
 bass = zeros(0,1);
 kass = zeros(0,1);
 offern = offer;
 return        % return
end
% Assign base and channel to new calls.
% Add some noise so the selection of base will be randomized.
mask = isnan(k)|isnan(b); % Assign only links without a channel.
bass = maxind(g+homargin/2*rand(size(g)),[],2);
bass(~mask) = nan; % Make all ongoing calls to nan.
[kass, offern] = radchanalloc(bass, offer);
reset = isnan(kass)&~mask; % Do not alter values that already was assigned.
bass(reset)=b(reset);
kass(reset)=k(reset);
bass = bass(:); % Make sure that input is only in first dimension.
kass = kass(:);
