function [ totald2dtocuInterference, d2dCapacityArray, avgD2DCapacity, d2dPowerAllocationMatrix, avgD2DCapacityConvergeProcess ] = mainAlgorithm( nochannels, noD2Dusers, d2dttod2drChannelGainMatrix, d2dttobaseChannelGainMatrix,...
    powerInitMatrix, piArray, cutod2drInterference, cuchannelAssignMatrix, P_max, Psub_max, noise )

% D2D-SCALE
%   Detailed explanation goes here



nOuter = 15; % number of times to update alpha and beta arrays
nInner = 150; % number of times to update lambda

%powerInitMatrix = 1/nochannels.*ones(nochannels,noD2Dusers);

lambda = zeros(nInner, noD2Dusers);
stepSize = 20;
[currentAlphaMatrix, currentBetaMatrix] = alphaBetaMatrixGen( nochannels, noD2Dusers, d2dttod2drChannelGainMatrix,...
    cutod2drInterference, cuchannelAssignMatrix, powerInitMatrix, noise );
alphaMatrix(:,:,1) = currentAlphaMatrix;
betaMatrix(:,:,1) = currentBetaMatrix;

lambda(1,:) = 1/(P_max*log(2)).*sum(currentAlphaMatrix)-min(piArray(ones(noD2Dusers,1),:)'.*d2dttobaseChannelGainMatrix);

for iOuter = 2:nOuter       
    for iInner = 1:nInner
        [ powerAllocateMatrix, sinrMatrix ] = fixedPoint( nochannels, noD2Dusers, d2dttod2drChannelGainMatrix, d2dttobaseChannelGainMatrix,...
            currentAlphaMatrix, lambda(iInner,:), piArray, cutod2drInterference, cuchannelAssignMatrix, noise);
        lambda(iInner+1,:) = max(0,lambda(iInner,:) + (stepSize/iInner).*(sum(powerAllocateMatrix) - P_max));
    end
    avgD2DCapacityConvergeProcess(iOuter - 1) = mean(sum(log2(1 + sinrMatrix)));
    [alphaMatrix(:,:,iOuter), betaMatrix(:,:,iOuter)] = alphaBetaMatrixGen( nochannels, noD2Dusers, d2dttod2drChannelGainMatrix,...
        cutod2drInterference, cuchannelAssignMatrix, powerAllocateMatrix, noise );
    currentAlphaMatrix = alphaMatrix(:,:,iOuter);
end

d2dCapacityArray = sum(log2(1 + sinrMatrix));
avgD2DCapacity = mean(d2dCapacityArray);
d2dPowerAllocationMatrix = powerAllocateMatrix;
totald2dtocuInterference = sum(d2dttobaseChannelGainMatrix.*d2dPowerAllocationMatrix,2);