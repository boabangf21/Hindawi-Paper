function [bho, kho, offern] = handoff(b, k, g, offer, homargin)
% DESCRIPTION [bho, kho, offern] = handoff(b, k, g, offer, homargin)
%  No handoff will be performed on b and k with values NaN. 
%  When base is selected an advantage is set to currently used
%  channel by homargin [dB].
% INPUT
%  b --         Base index for each link. Size is #links by one.
%               Elements of b that is equal to NaN are 
%               not going to be handed off. 
%  k --         Channel index for each link. Size is #links by one.
%               The value NaN indicates no hand off to this base. 
%  g --         Gain matrix with gain values in dB. Size is #links by #bases
%  offer --     Binary matrix. Zero indicate occupied channel. 
%               Size #bases by #channels.
%  homargin --  Standard deviation of the noise added to the
%               gain matrix when base selection is made. [dB]
% OUTPUT
%  bho --       As b but possibly with new base stations.
%  kass --      As k but possibly with new channels when
%               handoff was possible. Otherwise the old base and
%               channel would be kept.
%  offern --    Same size as offer but with new base station-
%               channel combination marked as occupied (zero)
%               and the old one released (set to one).
% TRY 
%  [b,k,offer]=handoff(1,1,[-100 -90],[0;1],0)
%   One link connected to base 1. Two bases, each with one 
%   channel are available. No noise is added.
%  [b,k,offer]=handoff(1,1,[-100 -90],[0 1;1 1],3)
%   One link connected to base 1. Two bases, each with two 
%   channels are available. 3 dB noise is added.

% by Magnus Almgren 000505

if isempty(b)   % the empty case
 bho = zeros(0,1);
 kho = zeros(0,1);
 offern = offer;
 break 
end

mask = ~isnan(k)&~isnan(b); % Do only act on already assigned channels.
m = (1:size(g,1))'; % make an index to each row in the gain matrix
ind1 = index(m(~isnan(b))',b(~isnan(b))',size(g));
g(ind1) = g(ind1) + homargin; % Give selected bases an advantage.
bhot = maxind(g,[],2); % desired bases
nohoff = (bhot==b)|isnan(b)|isnan(k); % No handoff on these. Keep as it is.
[khot, offern] = radchanalloc(bhot, offer); % kho is NaN where blocking occurred. 

% set as before
bho = b; 
kho = k; 

% Set new base and channel where hand off was successful.
hosucc = ~nohoff & b~=bhot & ~isnan(khot);
bho(hosucc) = bhot(hosucc); 
kho(hosucc) = khot(hosucc); 

% Update the allocation matrix.
offern = offer;
offern(index(  b(hosucc),  k(hosucc),size(offern))) = 1; % Set free.
offern(index(bho(hosucc),kho(hosucc),size(offern))) = 0; % Set occupied.

bho = bho(:); % Only spanning first dimension.
kho = kho(:);
