function obk=crechanplan(nb, nk, ncluster);
% DESCRIPTION obk=crechanplan(nb,nk,ncluster)
%  Makes a "frequency" plan with cluster size nb/ncluster.^2. 
%  The cells are numbered in the same way as in crecells. 
%  To get a correct symmetric frequency plan (nk*ncluster.^2)/nb 
%  has to evaluate into an integer.
% INPUT 
%  nb --        Number of base stations.
%  nk --        Number of channels in the system.
%  ncluster --  Square root of clusters.
% OUTPUT
%  obk --       Binary matrix, bases by channels.
%               obk(j,k)=1, if base j is an allocated channel.
% TRY          
%  spy(crechanplan(12,24,2))
% SEE ALSO     
%  crecells

nclu2 = ncluster.^2; % number of clusters
cellpclust = ceil(nb/nclu2); % clustersize;
obk = repmat(eye(cellpclust),[nclu2, ceil(nk/cellpclust)]);
obk = obk(1:nb,1:nk); % Correct the size if
                      % (nk*ncluster.^2)/nb was not an integer.



