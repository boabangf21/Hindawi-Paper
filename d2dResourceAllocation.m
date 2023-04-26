function [ totald2dtocuInterference, d2dCapacityArray, d2dPowerAllocationMatrix, powerConvergenceProcessEachChannel ] =d2dResourceAllocation( ...
    d2dttobaseChannelGainMatrix, d2dttod2drChannelGainMatrix, cutod2drInterference_rui15, ...
    P_max, Psub_max, noD2Dusers, nochannels, noise, piArray, cuchannelAssignMatrix )

nOuter = 500;
% initialize the transmit power price mu and the step size to update it
mu = zeros(1,noD2Dusers);
%mu(1,1) = 7.5;
stepSize = 0.03;
% number of iteration used to update the transmit power and interference
% price on each subchannel for each D2D pair
%powerConvergenceProcess = zeros(nochannels, noD2Dusers, maxIte);
%interferencePriceConvergenceProcess = zeros(nochannels, noD2Dusers, maxIte);
% initialize the transmit power
%powerInit = P_max/nochannels.*ones(1, noD2Dusers);
%powerCurrentArray = powerInit;
d2dCapacityEach = zeros(nochannels, noD2Dusers);
d2dPowerAllocationMatrix = zeros(nochannels, noD2Dusers);
% Loop from here
for iOuter = 1:nOuter %this is the loop for updating lagrangian multiplier
for iChannel = 1: nochannels %loop for each subchannel
% create the channel matrix for D2DT to D2DR for each iChannel, row
% is D2DT index, column is D2DR index
hArray = d2dttod2drChannelGainMatrix(:,:,iChannel);
diagHArray = reshape(diag(hArray),1,noD2Dusers);
% interference from cu on each subchannel, column is CU index, row is
% D2DR index
eachcutod2drInterferenceMatrix = cutod2drInterference_rui15(:,:,iChannel);
% reshape the interference from CU to D2DR, it is a new row of
% interference from the assigned CU to all D2DR on the iChannel
eachcutod2drInterference = eachcutod2drInterferenceMatrix(:,find(cuchannelAssignMatrix(iChannel,:)==1)).';
% stopping criterion
stopEps = 0.0001;
n0 = noise + eachcutod2drInterference;
% initialize power and max price initialization
%powerInit = stopEps.*ones(1, noD2Dusers);
powerInit=( P_max./nochannels).*ones(1, noD2Dusers);
powerCurrentArray = powerInit;
%currentSINR = (diagChannelMatrix.*powerCurrentArray)./(noise + eachcutod2drInterference + (powerCurrentArray*(channelMatrix)- diagChannelMatrix.*powerCurrentArray));
%priceCurrnetArray =(currentSINR.*currentSINR)./(powerCurrentArray.*diagChannelMatrix.*(1 + currentSINR).*log(2));
%priceCurrentArray = (diagChannelMatrix.*powerCurrentArray)./((1 + currentSINR).*log(2).*((noise + eachcutod2drInterference +(powerCurrentArray*(channelMatrix)- diagChannelMatrix.*powerCurrentArray)).^2));
%priceCurrentArray = ones(1, noD2Dusers)./n0;
iUpdate = 0;
powerConvergenceProcessiChannel = zeros(1, noD2Dusers);
%priceConvergenceProcessiChannel = zeros(1, noD2Dusers);
% maximal gap between two iterations
iteGapMax = 2*stopEps;
% get the piArray for each D2D pair
piIchannel = piArray(iChannel);
piArrayReshape = piIchannel(:,ones(1,noD2Dusers));
while (iteGapMax > stopEps) | (iUpdate < 10); %loop to update price and transmit power for each subchannel untill convergence
iUpdate = iUpdate + 1;
powerConvergenceProcessiChannel(iUpdate,:) = powerCurrentArray;
%priceConvergenceProcessiChannel(iUpdate,:) = priceCurrentArray;
% update the interference price according Eq. (10)
%priceCurrnetArray = (currentSINR).^2/(powerCurrentArray.*diagChannelMatrix.*(1 + currentSINR).*log(2));
% update the transmit power according equation (9), it is a row
% which each value is Eq. (9) for each D2D
equation9 = (1/log(2))./(piArrayReshape.*d2dttobaseChannelGainMatrix(iChannel,:) + mu(iOuter,:) )-( powerCurrentArray*(hArray)-powerCurrentArray.*diagHArray)./diagHArray; 
%equation9 = (1/log(2))./(priceCurrentArray*(hArray') - priceCurrentArray.*diagHArray) -(n0 + powerCurrentArray*(hArray)-powerCurrentArray.*diagHArray)./diagHArray;
powerArray(iUpdate, :) = max(min(equation9,Psub_max),0);
iteGapMax = max(abs((powerArray(iUpdate,:)- powerCurrentArray))./powerCurrentArray);
% update the currnet transmit power array
powerCurrentArray = powerArray(iUpdate, :);
currentSINR(iUpdate,:) = (powerCurrentArray.*diagHArray)./(n0 + powerCurrentArray*hArray-powerCurrentArray.*diagHArray);
% update the price on iChannel, it is row vector which each
% value is the price of each D2D
%priceArray(iUpdate,:) = 1./(log(2).*(n0 + powerCurrentArray*hArray-powerCurrentArray.*diagHArray));
%priceCurrentArray = priceArray(iUpdate,:);
%powerConvergenceProcess(iChannel,:,iUpdate) = reshape(powerCurrentArray,[1 nochannels 1]);
%interferencePriceConvergenceProcess(iChannel,:,iUpdate) = reshape(priceCurrnetArray,[1 nochannels 1]);
end
d2dCapacityEach(iChannel, :) = log2(1 + currentSINR(iUpdate, :));
d2dPowerAllocationMatrix(iChannel,:) = powerArray(iUpdate,:);
powerConvergenceProcessEachChannel{iChannel} = powerConvergenceProcessiChannel;
d2dt1PowerConvergenceProcess{iChannel} = powerConvergenceProcessiChannel(:,1);
end
mu(iOuter+1,:) = max(0, mu(iOuter,:) + (stepSize).*(sum(d2dPowerAllocationMatrix)-P_max));
end

% calculate the capacity of each D2D pair
d2dCapacityArray = sum(d2dCapacityEach);
% calculate the total interference form D2DT to cu on each subchannel, it
% is column vector, column is channel index
totald2dtocuInterference= sum(d2dttobaseChannelGainMatrix.*d2dPowerAllocationMatrix,2);
end