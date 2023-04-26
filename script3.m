

clear all;
nochannels = 16;
noD2Dusers = 2;
noCUusers = 2;

Dist = [100 200 300 400 500];%meters % increasing the  maximum distance between the D2D transmitters in order to  accommodate more D2D pairs and increase cell raduis

P_max =0.25;
%Pd2d_max = 0.4;
Psub_max =0.25;
%Pd2dsub_max = 0.4;
%noise power
noise = 0.1e-6;

%Pd2d_max = 0.1;
%Pd2dsub_max = 0.1;

% initialize the price for reusing each subchannel
%piArray = zeros(1,nochannels);
%piArray_rui15 = 0.5.*ones(1,nochannels);
% initialize the minimum data rate of each CU
%Rcu_min = 8;
%Rd2d_min = 70;
% initialize the number of realization
Nrealization =1;
deltaPi = 200000;
deltaPi_rui15 =200000;
transReciDisMin = 15;
transReciDisMax = 30;

%for iRealized = 1:nRealization

for id2d = 1:length(Dist)
    distance = Dist(id2d);
for iRealize = 1:Nrealization
        piArray = zeros(1,nochannels);
        piArray_rui15 = zeros(1,nochannels);
        piArray_centralized=zeros(1,nochannels);
        totald2dtocuInterference = zeros(nochannels, 1);
        totald2dtocuInterference_rui15 = zeros(nochannels, 1);
        totald2dtocuInterference_centralized=zeros(nochannels, 1);
        % Create the topology and channel gain matrix
        [ cutobaseChannelGainMatrix, d2dttobaseChannelGainMatrix, cutod2drChannelGainMatrix, ...
            d2dttod2drChannelGainMatrix] = channelGen3(noCUusers, noD2Dusers,distance, nochannels,transReciDisMin,transReciDisMax );

        % Assign channel to CU users
        [ newCUtoBaseChannelMatrix, cuchannelAssignMatrix ] = cuchannelAssign( nochannels, noCUusers, cutobaseChannelGainMatrix );

        % calculate the starting power allocation of CU
        [ cuCapacityArray, cuPowerAllocatedMatrix, cutod2drInterference, Target_interference ,Compare_newd2dInterference  ] = cuwaterfilling( P_max, Psub_max, cutobaseChannelGainMatrix, noise, ...
            cuchannelAssignMatrix, noCUusers, totald2dtocuInterference, nochannels, noD2Dusers, cutod2drChannelGainMatrix);
        cuCapacityMax = cuCapacityArray;
        cuCapacityConvergeProcess(1,:) = cuCapacityArray;
        
        
        [ cuCapacityArray_rui15, cuPowerAllocatedMatrix_rui15, cutod2drInterference_rui15,Target_interference_rui15,Compare_newd2dInterference_rui15 ] = cuwaterfilling_rui15( P_max, Psub_max, cutobaseChannelGainMatrix, noise, ...
            cuchannelAssignMatrix, noCUusers, totald2dtocuInterference_rui15, nochannels, noD2Dusers, cutod2drChannelGainMatrix);
       

        % calculate capacity of D2D by cellular mode
      %[ cellularModeCapacity(iRealize,id2d), d2dCellularModeCapacity, avgD2DCellularModeCapacity(iRealize,id2d) ] = cellularmode( P_max, Psub_max, cutobaseChannelGainMatrix, d2dttobaseChannelGainMatrix, noise,...
        %   noCUusers, nochannels, noD2Dusers);

        % Calculate capacity of D2D by heuristic algorithm
        [d2dCapacity_heuristic(iRealize,id2d), d2dPowerAllocationMatrix_heuristic, avgD2Dcapacity_heuristic(iRealize,id2d) ] =...
            heuristic(d2dttod2drChannelGainMatrix, cutod2drInterference,...
               P_max, Psub_max, noD2Dusers, nochannels, noise );
    
        % calculate the starting power allocation of D2DT
        [ totald2dtocuInterference, d2dCapacityArray, d2dPowerAllocationMatrix, powerConvergenceProcessEachChannel ] = d2dResourceAllocation( ...
            d2dttobaseChannelGainMatrix, d2dttod2drChannelGainMatrix, cutod2drInterference, ...
                P_max, Psub_max, noD2Dusers, nochannels, noise, piArray, cuchannelAssignMatrix );
                    Rd2d_min_rui15 = min(1/(2)*d2dCapacityArray);
                    
    %[totald2dtocuInterference_centralized, d2dCapacityArray_centralized, d2dPowerAllocationMatrix_centralized, powerConvergenceProcessEachChannel_centralized] = centralized( ...
   % d2dttobaseChannelGainMatrix, d2dttod2drChannelGainMatrix, cutod2drInterference_rui15, ...
 % P_max, Psub_max, noD2Dusers, nochannels,piArray_centralized, noise,  cuchannelAssignMatrix );
                    
                    
                    
        [ totald2dtocuInterference_rui15, d2dCapacityArray_rui15, d2dPowerAllocationMatrix_rui15, powerConvergenceProcessEachChannel_rui15 ] = rui15( ...
            d2dttobaseChannelGainMatrix, d2dttod2drChannelGainMatrix, cutod2drInterference, ...
                P_max, Psub_max, noD2Dusers, nochannels, noise, piArray_rui15, cuchannelAssignMatrix );
    
      
        %calculate again the capacity of CU users when there are interference from D2DTs        
        [ cuCapacityArray, cuPowerAllocatedMatrix, cutod2drInterference, Target_interference ,Compare_newd2dInterference ] = cuwaterfilling( P_max, Psub_max, cutobaseChannelGainMatrix, noise, ...
            cuchannelAssignMatrix, noCUusers, totald2dtocuInterference, nochannels, noD2Dusers, cutod2drChannelGainMatrix);
        cuCapacityMin = cuCapacityArray;
        cuCapacityConvergeProcess(2,:) = cuCapacityArray;
        Rcu_min =  min(cuCapacityMin + (cuCapacityMax-cuCapacityMin)./3);
   % Rcu_min = 80;
    
        [ cuCapacityArray_rui15, cuPowerAllocatedMatrix_rui15, cutod2drInterference_rui15,Target_interference_rui15,Compare_newd2dInterference_rui15 ] = cuwaterfilling_rui15( P_max, Psub_max, cutobaseChannelGainMatrix, noise, ...
            cuchannelAssignMatrix, noCUusers, totald2dtocuInterference_rui15, nochannels, noD2Dusers, cutod2drChannelGainMatrix);
        cuCapacityMin_rui15 = cuCapacityArray_rui15;
        Rcu_min_rui15 = min(cuCapacityMin_rui15 + (cuCapacityMax-cuCapacityMin_rui15)./3);
    %Rcu_min_rui15 =80;
        iPi = 0;
while any(cuCapacityArray < Rcu_min)
        %while any(Compare_newd2dInterference>mean( Target_interference))
        %%pricing interference scheme
            iPi = iPi + 1;
            %piArrayTemp = piArray(iPi - 1);
            cuindex = cuCapacityArray < Rcu_min;
            %cuindex = Compare_newd2dInterference> mean( Target_interference);
            cuchannelAssignMatrix(:,cuindex);
            if sum(cuindex) == 1
                deltaPiArray = deltaPi.*cuchannelAssignMatrix(:,cuindex)';
            else deltaPiArray = sum(deltaPi.*cuchannelAssignMatrix(:,cuindex)');
            end
            piArray = piArray + deltaPiArray;
        
            [ totald2dtocuInterference, d2dCapacityArray, d2dPowerAllocationMatrix, powerConvergenceProcessEachChannel ] = d2dResourceAllocation( ...
                d2dttobaseChannelGainMatrix, d2dttod2drChannelGainMatrix, cutod2drInterference, ...
                    P_max, Psub_max, noD2Dusers, nochannels, noise, piArray, cuchannelAssignMatrix );



            [ cuCapacityArray, cuPowerAllocatedMatrix, cutod2drInterference,Target_interference ,Compare_newd2dInterference ] = cuwaterfilling( P_max, Psub_max, cutobaseChannelGainMatrix, noise, ...
                cuchannelAssignMatrix, noCUusers, totald2dtocuInterference, nochannels, noD2Dusers, cutod2drChannelGainMatrix);
            cuCapacityConvergeProcess(iPi+1,:) = cuCapacityArray;
            interferenceConvergeProcess(iPi+1,:)=Compare_newd2dInterference;
         %d2dCapacityConverge(iPi+1,:)=d2dCapacityArray;
end
        avgD2DCapacity(iRealize,id2d) = mean(d2dCapacityArray)

 
                  
  iPd=0;
    
      while any(cuCapacityArray_rui15 < Rcu_min_rui15)
            iPd=iPd+1;
           cuindex = cuCapacityArray_rui15 < Rcu_min_rui15;
            %while any(Compare_newd2dInterference>mean( Target_interference))
                  %cuindex = Compare_newd2dInterference> mean( Target_interference);
            if sum(cuindex) == 1
                deltaPiArray_rui15 = deltaPi_rui15.*cuchannelAssignMatrix(:,cuindex)';
            else deltaPiArray_rui15 = sum(deltaPi_rui15.*cuchannelAssignMatrix(:,cuindex)');
            end
            piArray_rui15 = piArray_rui15 + deltaPiArray_rui15;
        
            [ totald2dtocuInterference_rui15, d2dCapacityArray_rui15, d2dPowerAllocationMatrix_rui15, powerConvergenceProcessEachChannel_rui15 ] = rui15( ...
                d2dttobaseChannelGainMatrix, d2dttod2drChannelGainMatrix, cutod2drInterference_rui15,  ...
                    P_max, Psub_max, noD2Dusers, nochannels, noise, piArray_rui15, cuchannelAssignMatrix);
                
        
            [ cuCapacityArray_rui15, cuPowerAllocatedMatrix_rui15, cutod2drInterference_rui15,Target_interference_rui15,Compare_newd2dInterference_rui15 ] = cuwaterfilling_rui15( P_max, Psub_max, cutobaseChannelGainMatrix, noise, ...
                cuchannelAssignMatrix, noCUusers, totald2dtocuInterference_rui15, nochannels, noD2Dusers, cutod2drChannelGainMatrix);
        interferenceConvergeProcess_rui15(iPd+1,:)= Compare_newd2dInterference_rui15;
            cuCapacityConvergeProcess_rui15(iPd+1,:) = cuCapacityArray_rui15;
    %d2dCapacityConverge_rui15(iPd+1,:)= d2dCapacityArray_rui15;
      end
        avgD2DCapacity_rui15(iRealize,id2d) = mean(d2dCapacityArray_rui15)
 end
end


