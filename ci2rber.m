function rber = ci2rber(ci)
% DESCRIPTION res = ci2rber(x)
% a simple model for Raw Bit Error Rate
% INPUT 
% ci --   Carrier to Interference in dB
% OUTPUT 
% rber -- Raw Bit Error Rate
% TRY 
% plot(-30:10, ci2rber(-30:10))
rber = (1 + erf(-db2lin(ci/2)))/2;
