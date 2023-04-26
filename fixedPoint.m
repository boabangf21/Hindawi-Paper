function [ powerAllocateMatrix, sinrMatrix ] = fixedPoint( nochannels, noD2Dusers, d2dttod2drChannelGainMatrix, d2dttobaseChannelGainMatrix,...
    alphaMatrix, lambda, piArray, cutod2drInterference, cuchannelAssignMatrix, noise)
% This function return the powerAllocateMatrix of d2dt for each subchannel
% as the solution of Eq. 12, row is subchannel index, column is d2dt index
%   sinrMatrix: sinr matrix, row is subchannel, column is d2dt index

for iChannel = 1:nochannels
    % get the channel matrix of d2dt->d2dr for each subchannel
    channelArray = d2dttod2drChannelGainMatrix(:,:,iChannel);
    diagChannelArray = reshape(diag(channelArray),1,noD2Dusers);
    % get the channel array of d2dt->base for each subchannel
    d2dttobaseChannelGainArray = d2dttobaseChannelGainMatrix(iChannel,:);
    alphaArray = alphaMatrix(iChannel,:);
    powerInit = 0.5.*ones(1,noD2Dusers);
    % get the interference from cu->d2dr for each subchannel
    eachcutod2drInterferenceMatrix = cutod2drInterference(:,:,iChannel);
    eachcutod2drInterference = eachcutod2drInterferenceMatrix(:,find(cuchannelAssignMatrix(iChannel,:)==1)).';
    
    piIchannel = piArray(iChannel);
    piArrayReshape = piIchannel(:,ones(1,noD2Dusers));
        
    powerCurrentArray = powerInit;
    powerConvergenceProcessArray(1,:) = powerCurrentArray;
    for i = 2:10
        powerConvergenceProcessArray(i,:) = interferenceFunction( channelArray,...
            d2dttobaseChannelGainArray, alphaArray, powerCurrentArray, noD2Dusers, lambda, piArrayReshape, eachcutod2drInterference, noise );
        powerCurrentArray = powerConvergenceProcessArray(i,:);
    end
    powerAllocateMatrix(iChannel,:) = powerConvergenceProcessArray(i,:);
    sinrMatrix(iChannel,:) = powerCurrentArray.*diagChannelArray./(noise + powerCurrentArray*channelArray - powerCurrentArray.*diagChannelArray + eachcutod2drInterference);
end


