function [ interferenceFunctionValue ] = interferenceFunction( channelArray,...
    d2dttobaseChannelGainArray, alphaArray, d2dTransmitPowerArray, noD2Dusers, lambda, piArray, eachcutod2drInterference,noise )
%Interference function value for each subchannel, it is a row vector

%   hMatrix: Channel gain matrix d2dt->d2dr for each subchannel
%   d2dttobaseChannelGainArray: row vector of channel gains from d2dt->base
%   alphaArray: row vector of alpha for each subchannel
%   d2dTransmitPowerMatrix: row vector of d2dt transmit power
%   cutod2drInterference: row vector of interference from cu to d2dr for
%   each subchannel

% create the row channel gain between d2dt->d2dr of a d2d pair
diagChannelArray = reshape(diag(channelArray),1,noD2Dusers);

iplusn = noise + d2dTransmitPowerArray*channelArray - d2dTransmitPowerArray.*diagChannelArray + eachcutod2drInterference;

%inr = noise + d2dTransmitPowerArray*hMatrix - d2dTransmitPowerArray.*diaghMatrix + eachcutod2drInterference;

% calculate alpha^n_i*h^n_di in Eq. (12)
alphaHMatrix = alphaArray(ones(noD2Dusers,1),:).*channelArray;

A = alphaHMatrix./iplusn(ones(noD2Dusers,1),:);
B = sum(A,2)' - reshape(diag(A),1,noD2Dusers);

% calculate equation 12
eq12 = alphaArray./((piArray.*d2dttobaseChannelGainArray + lambda).*log(2) + B);

interferenceFunctionValue = eq12;







