function nmob = mrequest(sofftraf, slope)
% DESCRIPTION nmob = mrequest(sofftraf, slope)
%  Generate new users requesting a base and channel assignment.
%  veclength is the number of inactive users from which we
%  select nmobavg users in average that are going to be active.
%  100 is the minimum pool of inactive users when nmobavg is small.
% INPUT 
%  sofftraf -- The traffic offered to the simulated system.
%  slope --    Ratio dt to mht.
% OUTPUT
%  nmob --     A realisation of number of new users to create.
% TRY      
%  mrequest(1000, 0.01)
%  mrequest(1000, inf)
% SEE ALSO 
%  terminate

% by Magnus Almgren 000508

nmobavg = sofftraf *(1-exp(-slope)); 
veclength = round(max(100, 10*nmobavg)); 
nmob = sum(rand(veclength,1) < nmobavg/veclength);

