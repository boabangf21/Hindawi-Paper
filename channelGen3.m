function [ cutobaseChannelGainMatrix, d2dttobaseChannelGainMatrix, cutod2drChannelGainMatrix, d2dttod2drChannelGainMatrix] = channelGen3(noCUusers, noD2Dusers,distance, nochannels,transReciDisMin,transReciDisMax )

%=========================================================================
%   This function returns the channel gain matrix in linear
%=========================================================================
% cutobaseChannelGainMatrix is channel gain matrix from CU to BS, row is
% channel index, column is CU index

% d2dttobaseChannelGainMatrix is channel gain between D2DT and BS, row is
% channel index, column is D2DT index

% cutod2drChannelGainMatrix is channel gain between cu to d2dr, row is
% D2DR, column is CU and page is channel index

% d2dttod2drChannelGainMatrix is channel gain between D2DT and D2DR, row is
% D2DT, column is D2DR and page is channel index


% Create cells, xyb is the coordination of the BS
%nochannels = 16;
%noCUusers = 2;
%noD2Dusers = 200;
cellRadiusMax = 500; cellRadiusMin = 200;
sps = 1; km = 1; lm = 0; ncluster = 1;
%noOfBaseStation = km^2 + lm^2 + km*lm;
%logNormalMean = 0; logNormalDeviation = 8;
[xyb, fib, rombvec] = crecells(cellRadiusMax, sps, km, lm, ncluster);
temp = xyb;
noOfBaseStation = km^2 + lm^2 + km*lm;
%-------------------------------------------------------------------------

% create the location of D2DT and D2DR
% Minimum range between D2DT and D2DR is set to 15m
% Maximum range between D2DT and D2DR is set to 30m
%transReciDisMin = 15;
%transReciDisMax = 30;

% create 2000 random locations for D2DT
disMobile2Base = cellRadiusMin + (cellRadiusMax - cellRadiusMin)*rand(1,2000);
tempphi = rand(1,2000)*2*pi;
mobileCoordinateRow = complex(disMobile2Base.*cos(tempphi), disMobile2Base.*sin(tempphi));
mobileCoordinateColumn = mobileCoordinateRow.';
mobileDistancePair = abs(mobileCoordinateRow(ones(2000,1),:)-mobileCoordinateColumn(:,ones(1,2000)));
% create random location of CUs
index = randperm(2000);
cuIndex = index(1:noCUusers);
xycu = mobileCoordinateColumn(cuIndex,1);
xycurow = xycu.';

% choose location of D2DT candidate which distance to CU is larger than
% 100m and satisfies distance condition to BS
cuandd2dDisCandidate = mobileDistancePair(cuIndex,:);
d2dtCandidateIndex = find(sum(cuandd2dDisCandidate >=100) == noCUusers); % this should be work with more than one CU user
% Obtain the coordinate of candidate D2DTs
d2dtCandidateCoordinateRow = mobileCoordinateRow(d2dtCandidateIndex);
d2dtCandidateCoordinateColumn = d2dtCandidateCoordinateRow.';
% randomly choose the locaion for the first D2DT
nod2dtCandidate = length(d2dtCandidateCoordinateRow);
newIndex = randperm(length(d2dtCandidateCoordinateRow));
xyd2dt = zeros(noD2Dusers,1);
xyd2dt(1) = d2dtCandidateCoordinateColumn(newIndex(1));
% randomly choose the subsequent location of other D2DTs which the
% distances between any 2 D2DTs is larger than 100
for i = 2:noD2Dusers
    previousD2DSet = xyd2dt(1:i-1,:);
    tempd2dtod2dDisPair = abs(previousD2DSet(:,ones(1,nod2dtCandidate))-d2dtCandidateCoordinateRow(ones(i-1,1),:));
    if i == 2
        tempd2dCandidateIndex =find((tempd2dtod2dDisPair<=distance));
        %tempd2dCandidateIndex = find((tempd2dtod2dDisPair >= 50)&(tempd2dtod2dDisPair <= 2000));
    else
         tempd2dCandidateIndex =find(sum((tempd2dtod2dDisPair<=distance))== i-1);
        %tempd2dCandidateIndex = find(sum((tempd2dtod2dDisPair >= 50)&(tempd2dtod2dDisPair <= 2000)) == i-1);
    end
    tempd2dIndex = randperm(length(tempd2dCandidateIndex));
    xyd2dt(i) = d2dtCandidateCoordinateColumn(tempd2dCandidateIndex(tempd2dIndex(1)));
end

% generate the location of D2DR randomly
for j = 1:noD2Dusers
    distanceToD2DR = transReciDisMin + (transReciDisMax - transReciDisMin)*rand(1,1);
    phi = rand(1,1)*2*pi;
    coordinateD2DRtoBase = xyd2dt(j) + complex(distanceToD2DR*cos(phi), distanceToD2DR*sin(phi));
    if mod(j,noOfBaseStation) ~= 0
        xyd2dr(j) = temp(mod(j,noOfBaseStation)) + coordinateD2DRtoBase;
    else
        xyd2dr(j) = temp(noOfBaseStation) + coordinateD2DRtoBase;
    end
end
xyd2dr = xyd2dr.';

% create the distance array between CUs and BS
cutobaseDistanceArray = abs(xycu-xyb).';

% create the distance matrix between D2DTs and D2DRs
xyd2drRow = xyd2dr.';
d2dtod2drDistanceMatrix = abs(xyd2dt(:,ones(1,noD2Dusers)) - xyd2drRow(ones(1,noD2Dusers),:));

% create the distance matrix between CUs and D2DRs, row is CU index, column
% is D2DR index
cutod2drDistanceMatrix = abs(xycurow(ones(noD2Dusers,1),:) - xyd2dr(:,ones(1,noCUusers)));

% create the distance array between D2DTs and BS
d2dttobaseDistanceArray = abs(xyd2dt - xyb).';

% create the channel gain array between CUs and BS in linear unit (not DB),
% row is subchannel index, column is CU index
% channel gain is normalized by noise which is 170dBm - 30 = 140 dbm
% db2lin change unit from dBm to normal
for i = 1:nochannels
    A = normrnd(0,4,25, noCUusers);
    shadowing1(i,:) = mean(A);
    B = normrnd(0,4,25, noD2Dusers);
    shadowing2(i,:) = mean(B);
    for j = 1:noD2Dusers
        C = normrnd(0,4,25, noCUusers);
        c = mean(C);
        shadowing3(j,:,i) = c;
        D = normrnd(0,4,25, noD2Dusers);
        shadowing4(j,:,i) = mean(D);
    end
end  
cutobaseChannelGainMatrix = db2lin(-37 - 30.*log10(cutobaseDistanceArray(ones(nochannels,1),:)) - shadowing1 + 65);

% create the channel gain matrix between D2DT to base station, row is
% subchannel index, column is D2DT index
%for i = 1:nochannels
%    B = normrnd(0,4,50, noD2Dusers);
%    shadowing2(i,:) = mean(B);
%end  
d2dttobaseChannelGainMatrix = db2lin(-37 - 30.*log10(d2dttobaseDistanceArray(ones(nochannels,1),:)) - shadowing2 + 65);

%create 3 dimensional channel gain matrix between CU to D2DR, row is D2DR
%index, column is CU index, page is subchannel index
%for i = 1:nochannels
 %   for j = 1:noD2Dusers
%        C = normrnd(0,4,50, noCUusers);
 %       c = mean(C);
 %       shadowing3(j,:,i) = c;
 %   end
%end
cutod2drChannelGainMatrix = db2lin(-37 - 30.*repmat(log10(cutod2drDistanceMatrix),[1 1 nochannels]) - shadowing3 + 65);

% create 3 dimensional channel gain matrix between D2DTs and D2DRs, row is
% D2DT index, column is D2DR index, page is subchannel index
%for i = 1:nochannels
%    for j = 1:noD2Dusers
%        D = normrnd(0,4,50, noD2Dusers);
%        shadowing4(j,:,i) = mean(D);
%    end
%end
d2dttod2drChannelGainMatrix = db2lin(-37 - 30.*repmat(log10(d2dtod2drDistanceMatrix),[1 1 nochannels]) - shadowing4 + 65);
% Caculate the real distance between BSs and users
% n is the number of channel realization for shadowing
% channelGain is channel gain matrix, represent pathloss and log-normal shadowing, row is BS index and column is user
% index
%n = 100;
%channelGain = zeros(noOfBaseStation, noCUusers);
%for i = 1:noOfBaseStation
%    distancePairTemp = reshape(abs(xyb(ones(noCUusers,1),i)-xym),1,noCUusers);
%    distancePair(i,:) = distancePairTemp;
%    shadowing = normrnd(logNormalMean, logNormalDeviation, n, 1);
%    channelGain(i,:) = -33 - 36.7*log10(distancePair(i,:)) + mean(shadowing);
%end

%xyv = zeros(noOfUsers,1);

%lobevector = zeros(360,1);
%corrdist = 110;
%attenconst = -21;
%[lognmap, mapvec] = crelognmap(xyb, rombvec, corrdist);
%raa = 0.5;
% Create the channel gain in dB
%g = pathgain(xym, xyb, fib, rombvec, attenconst, attenFactor, logNormalDeviation, raa, lobevector, lognmap, mapvec);

% Return the channel gain (normalized by the noise which is 170 dBm - 30 = 140 dB) in linear
%hArray = db2lin((channelGain + 130));
end

