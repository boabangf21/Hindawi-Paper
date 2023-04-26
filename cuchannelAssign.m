function [ newCUtoBaseChannelMatrix, channelAssignMatrix ] = cuchannelAssign( nochannels, noCUusers, cutobaseChannelGainMatrix )
%This function assign subchannel to CU users using a round-robin scheme

% cutobaseChannelGainMatrix is the channel gain matrix between CU to base
% (nochannels x noCUusers)

%   newCUtoBaseChannelMatrix is new channel gain matrix (nochannels x
%   noCUusers) in which only the channel gain of the assigned in each
%   column

% channelAssignMatrix is the channel assign matrix (nochannels x noCUusers)
% 1 means assigned, 0 means not assigned

channelAssignMatrixTemp = zeros(nochannels, noCUusers);
channelGainTemp         = cutobaseChannelGainMatrix;

partition = 1:noCUusers:nochannels;

for n = 1:length(partition)-1
    for k = 0:noCUusers - 1
        [a, index] = max(channelGainTemp(partition(n)+k,:));
        channelAssignMatrixTemp(partition(n)+k, index) = 1;
        channelGainTemp(partition(n)+k:partition(n)+noCUusers-1, index) = 0;
    end
end

for i = partition(n+1)
    for j = 0:(nochannels - i)
        [a, index] = max(channelGainTemp(i+j,:));
        channelAssignMatrixTemp(i+j, index) = 1;
        channelGainTemp(i+j:nochannels, index) = 0;
    end
end


%for n = 1:noCUusers:(nochannels - 1)
%    for k = 0:noCUusers-1
%        [a, index] = max(channelGainTemp(n+k,:));
%        channelAssignMatrixTemp(n+k, index) = 1;
%        channelGainTemp(n+k:n+noCUusers-1, index) = 0;
%    end
%end
%[a, index] = max(channelGainTemp(nochannels,:));
%channelAssignMatrixTemp(nochannels, index) = 1;

channelAssignMatrix         = channelAssignMatrixTemp;
newCUtoBaseChannelMatrix    = channelAssignMatrix.*cutobaseChannelGainMatrix;