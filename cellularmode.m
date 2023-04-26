function [ cellularModeCapacity, d2dCellularModeCapacity, avgD2DCellularModeCapacity ] = cellularmode( P_max, Psub_max, cutobaseChannelGainMatrix, d2dttobaseChannelGainMatrix, noise,...
    noCUusers, nochannels, noD2Dusers)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


cellularModeChannelGainMatrix = [cutobaseChannelGainMatrix d2dttobaseChannelGainMatrix];

% assign channel to CU users
cellularModenoCUusers = noCUusers + noD2Dusers;

channelAssignMatrixTemp = zeros(nochannels, noCUusers);
channelGainTemp = cellularModeChannelGainMatrix;
partition = 1:cellularModenoCUusers:nochannels;

for n = 1:length(partition)-1
    for k = 0:cellularModenoCUusers - 1
        [a, index] = max(channelGainTemp(partition(n)+k,:));
        channelAssignMatrixTemp(partition(n)+k, index) = 1;
        channelGainTemp(partition(n)+k:partition(n)+cellularModenoCUusers-1, index) = 0;
    end
end

for i = partition(n+1)
    for j = 0:(nochannels - i)
        [a, index] = max(channelGainTemp(i+j,:));
        channelAssignMatrixTemp(i+j, index) = 1;
        channelGainTemp(i+j:nochannels, index) = 0;
    end
end
cellularModeChannelAssignMatrix = channelAssignMatrixTemp;

% perform waterfilling algorithm with CU users
for n = 1:cellularModenoCUusers
    cuchannelAssignedLength = sum(cellularModeChannelAssignMatrix(:,n));
    index = find(cellularModeChannelAssignMatrix(:,n) == 1);
    newchannelArray         = cellularModeChannelGainMatrix(index,n);
    % calculate the column on interference from D2DT to cu on each subchannel 
    %newd2dInterference = totald2dtocuInterference(index, n);
    noisetosubchannel       = noise./newchannelArray;
        
    initPowerAllo = min(Psub_max, (P_max + sum(noisetosubchannel))./cuchannelAssignedLength - noisetosubchannel);
    
    while(length( find(initPowerAllo < 0 )) > 0 )
        negIndex       = find(initPowerAllo <= 0);
        posIndex       = find(initPowerAllo >  0);
        nSubchannelRem = length(posIndex);
        initPowerAllo(negIndex) = 0;
        
        snrRem                  = noisetosubchannel(posIndex);
        powerAlloTemp           = min(Psub_max, (P_max + sum(snrRem))./nSubchannelRem - snrRem);
        initPowerAllo(posIndex) = powerAlloTemp;
    end
    cellularModePowerAllocatedMatrix(index, n)    = initPowerAllo;
    cellularModeCapacityArray(n)          = sum(log2(1 + initPowerAllo./noisetosubchannel));
end
cellularModeCapacity = sum(1+cellularModeCapacityArray);
d2dCellularModeCapacity = sum(cellularModeCapacityArray(noCUusers+1:noCUusers+noD2Dusers));
avgD2DCellularModeCapacity = d2dCellularModeCapacity/noD2Dusers;
end

