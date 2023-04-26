function [ cuCapacityArray_rui15, cuPowerAllocatedMatrix_rui15, cutod2drInterference_rui15,Target_interference_rui15,Compare_newd2dInterference_rui15 ] = cuwaterfilling_rui15( P_max, Psub_max, channelGainMatrix, noise, ...
    cuchannelAssign, noCUusers, totald2dtocuInterference_rui15, nochannels, noD2Dusers, cutod2drChannelGainMatrix)
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


cutod2drInterference_rui15 = zeros(noD2Dusers, noCUusers, nochannels);
totald2dtocuInterference_rui15 = totald2dtocuInterference_rui15(:, ones(1,noCUusers));
for n = 1:noCUusers
    cuchannelAssignedLength = sum(cuchannelAssign(:,n));
    index = find(cuchannelAssign(:,n) == 1);
    newchannelArray         = channelGainMatrix(index,n);
    % calculate the column on interference from D2DT to cu on each subchannel 
    newd2dInterference_rui15 = totald2dtocuInterference_rui15(index, n) ;
   % Compare_newd2dInterference(n)=sum(newd2dInterference + noise);
    noisetosubchannel_3       = (noise + newd2dInterference_rui15)./newchannelArray;
    Compare_newd2dInterference_rui15(n)=sum(noisetosubchannel_3);
     noisetosubchannel_2   = (noise)./newchannelArray;   
    initPowerAllo = min(Psub_max, (P_max + sum(noisetosubchannel_2))./cuchannelAssignedLength - noisetosubchannel_2);
    %initPowerAllo = min(Psub_max)./(cuchannelAssignedLength);
    while(length( find(initPowerAllo < 0 )) > 0 )
        negIndex       = find(initPowerAllo <= 0);
        posIndex       = find(initPowerAllo >  0);
        nSubchannelRem = length(posIndex);
        initPowerAllo(negIndex) = 0;
        
        snrRem                  = noisetosubchannel_3(posIndex);
        powerAlloTemp           = min(Psub_max, (P_max + sum(snrRem))./nSubchannelRem - snrRem);
        initPowerAllo(posIndex) = powerAlloTemp;
    end
    cuPowerAllocatedMatrix_rui15(index, n)    = initPowerAllo;
    cuCapacityArray_rui15(n)          = sum(log2(1 + initPowerAllo./noisetosubchannel_3));
    % calculate the interference to D2D users
    B = initPowerAllo.';
    C = B(ones(noD2Dusers,1),:);
    A = cutod2drChannelGainMatrix(:,n,index).*reshape(C,[noD2Dusers 1 length(index)]);
    cutod2drInterference_rui15(:, n, index)           = A;
    
    Target_interference_rui15= (initPowerAllo)./((2.^( cuCapacityArray_rui15(n))-1)-noise);
end

