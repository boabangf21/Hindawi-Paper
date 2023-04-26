function [d2dCapacity_heuristic, d2dPowerAllocationMatrix_heuristic, avgD2Dcapacity_heuristic ] =...
    heuristic(d2dttod2drChannelGainMatrix, cutod2drInterference,...
    P_max, Psub_max, noD2Dusers, nochannels, noise )

%Heuristic algorithm in Zulhasnine10
%   Detailed explanation goes here

% create the channel gain matrix of d2d pair, row is channel index, column
% is d2d pair index
d2dpairChannelGainMatrix = zeros(nochannels, noD2Dusers);
for i = 1:nochannels
    d2dpairChannelGainMatrix(i,:) = reshape(diag(d2dttod2drChannelGainMatrix(:,:,i)),1,noD2Dusers);
end

% assign channel for d2d pair
[ newCUtoBaseChannelMatrix, d2dChannelAssignMatrix ] = cuchannelAssign( nochannels, noD2Dusers, d2dpairChannelGainMatrix );

for j = 1:nochannels
    A = sum(cutod2drInterference(:,:,j),2);
    cutod2drInterferenceMatrix(j,:) = d2dChannelAssignMatrix(j,:).*A.';
end
    

for n = 1:noD2Dusers
    d2dchannelAssignedLength = sum(d2dChannelAssignMatrix(:,n));
    index = find(d2dChannelAssignMatrix(:,n) == 1);
    newchannelArray         = d2dpairChannelGainMatrix(index,n);
    % calculate the column of interference from cu to d2dr on each
    % subchannel
    %cutod2drInterferenceMatrix = cutod2drInterference(:,n,index);
    newcutod2drInterference = cutod2drInterferenceMatrix(index,n);
    noisetosubchannel_1       = (noise + newcutod2drInterference)./newchannelArray;
    noisetosubchannel   = (noise )./newchannelArray;
        
    initPowerAllo = min(Psub_max, (P_max + sum(noisetosubchannel))./d2dchannelAssignedLength - noisetosubchannel);
    
    while(length( find(initPowerAllo < 0 )) > 0 )
        negIndex       = find(initPowerAllo <= 0);
        posIndex       = find(initPowerAllo >  0);
        nSubchannelRem = length(posIndex);
        initPowerAllo(negIndex) = 0;
        
        snrRem                  = noisetosubchannel(posIndex);
        powerAlloTemp           = min(Psub_max, (P_max + sum(snrRem))./nSubchannelRem - snrRem);
        initPowerAllo(posIndex) = powerAlloTemp;
    end
    heuristicPowerAllocatedMatrix(index, n)    = initPowerAllo;
    heuristicD2DCapacityArray(n)          = sum(log2(initPowerAllo./noisetosubchannel_1));
    
end
d2dCapacity_heuristic = sum(heuristicD2DCapacityArray);
avgD2Dcapacity_heuristic = d2dCapacity_heuristic./noD2Dusers;
d2dPowerAllocationMatrix_heuristic = heuristicPowerAllocatedMatrix;
end

