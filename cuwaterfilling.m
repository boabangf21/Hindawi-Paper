function [ cuCapacityArray, cuPowerAllocatedMatrix, cutod2drInterference,Target_interference,Compare_newd2dInterference ] = cuwaterfilling( P_max, Psub_max, channelGainMatrix, noise, ...
    cuchannelAssign, noCUusers, totald2dtocuInterference, nochannels, noD2Dusers, cutod2drChannelGainMatrix)
%Waterfilling algorithm for CUs
%   P_max is the total transmit power for each CU
%   channelGainMatrix is the channel gain matrix between CU to BS
%   (nochannels x noCUusers)

%   noise is the noise power
%   cuchannelAssign is the channel assignment (nochannels x noCUusers)
%   noCUusers is the number of CU users

%   d2dInterference is the total interference on each subchannel from D2DT to base station
%   (row vector)

%   nochannels is the number of subchannels
%   cuCapacityArray is a row represents capacity of each CU
%   cuPowerAllocatedMatrix is the matrix of CU power allocation, row is
%   subchannel index, column is CU index, 0 means that subchannel is not
%   assigned

%   cutod2drInterference is 3 dimensional matrix, column is CU index, row
%   is D2D index, page is subchannel index


cutod2drInterference = zeros(noD2Dusers, noCUusers, nochannels);
totald2dtocuInterference = totald2dtocuInterference(:, ones(1,noCUusers));
for n = 1:noCUusers
    cuchannelAssignedLength = sum(cuchannelAssign(:,n));
    index = find(cuchannelAssign(:,n) == 1);
    newchannelArray         = channelGainMatrix(index,n);
    % calculate the column on interference from D2DT to cu on each subchannel 
    newd2dInterference = totald2dtocuInterference(index, n) ;
   % Compare_newd2dInterference(n)=sum(newd2dInterference + noise);
    noisetosubchannel       = (noise + newd2dInterference)./newchannelArray;
    Compare_newd2dInterference(n)=sum(noisetosubchannel);
     noisetosubchannel_1   = (noise)./newchannelArray;   
    initPowerAllo = min(Psub_max, (P_max + sum(noisetosubchannel_1))./cuchannelAssignedLength - noisetosubchannel_1);
    %initPowerAllo = min(Psub_max)./(cuchannelAssignedLength);
    while(length( find(initPowerAllo < 0 )) > 0 )
        negIndex       = find(initPowerAllo <= 0);
        posIndex       = find(initPowerAllo >  0);
        nSubchannelRem = length(posIndex);
        initPowerAllo(negIndex) = 0;
        
        snrRem                  = noisetosubchannel(posIndex);
        powerAlloTemp           = min(Psub_max, (P_max + sum(snrRem))./nSubchannelRem - snrRem);
        initPowerAllo(posIndex) = powerAlloTemp;
    end
    cuPowerAllocatedMatrix(index, n)    = initPowerAllo;
    cuCapacityArray(n)          = sum(log2(1 + initPowerAllo./noisetosubchannel));
    % calculate the interference to D2D users
    B = initPowerAllo.';
    C = B(ones(noD2Dusers,1),:);
    A = cutod2drChannelGainMatrix(:,n,index).*reshape(C,[noD2Dusers 1 length(index)]);
    cutod2drInterference(:, n, index)           = A;
    
    Target_interference= (initPowerAllo)./((2.^( cuCapacityArray(n))-1)-noise);
end

