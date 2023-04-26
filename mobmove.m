function [xym,xyv] = mobmove(xym,xyv,amean,vmean,dt,rombvec)
% DESCRIPTION [xym,xyv] = mobmove(xym,xyv,amean,vmean,dt,rombvec)
%  Move and initiate position of users.
% INPUT
%  xym --     Complex position of users possibly. NaN for users that have not yet
%             received any position and thus are regarded as new.
%  xyv --     Complex speed [m/s].
%  amean --   Average acceleration of mobiles.
%  vmean --   Mean speed of mobiles. 
%  dt --      Time interval since last call to this routine.
%  rombvec -- Two complex vectors determining the area within which 
%             mobiles are folded.
% OUTPUT
%  xym --     New position. 
%  xyv --     New speed.
% TRY  
%  plot(mobmove(zeros(1000,1),zeros(1000,1),1,10,1,[1 i]*1000),'.')
%  plot(mobmove(zeros(1000,1),zeros(1000,1),1,10,1,[1 i]*1),'.')
%  plot(mobmove(nans(1000,1),zeros(1000,1),1,10,1,[1+i i]*1000),'.')
% SEE ALSO 
%  mrequest, terminate

% by Magnus Almgren 000505

% Move users.
raa = exp(-mdiv(dt.*amean,vmean));
% raa is the correlation between the new and the old xyv
% and does not impact the variance.
xyv = xyv.*raa+sqrt(1-raa.^2).*vmean/sqrt(pi/2).*irandn(size(xyv));
xym = xym + xyv.*dt;

% Give new mobiles speed and position.
newones = isnan(xym);

% Generate position and wrap into simulated area.
xym(newones)= 0;
if exist('rombvec','var')
 if all(isfinite(rombvec))
  xym(newones) = rand(size(xyv(newones),1),2)*rombvec(:);
  xym = wrapinto(xym+mean(rombvec),rombvec,'romb')-mean(rombvec);
 end
end
% Init speed.
xyv(newones) = vmean/sqrt(pi/2).*irandn(size(xyv(newones))); 

xym = xym(:); % Make sure only first dimension.
xyv = xyv(:);
