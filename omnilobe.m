function lobevector = omnilobe(gain)
% DESCRIPTION lobevector = omnilobe
%  Creates an antenna lobe vector. The antenna lobe is a
%  circle with radius 'gain' dB, if 'gain' is present, otherwise
%  with radius 3.0244 dB. The gain 3.0244 is chosen because
%  its sum in all direction will be equal to the sum in the
%  sinclobe case.
% OUTPUT
%  lobevector --  A row vector containing the gain in each 1 degree direction.
%                 The length of the vector is 360 where last index corresponds
%                 zero degrees and first index to 1 degree.
% TRY
%  cd ../examples; antenna_lobe; cd ../runemain
%  sum(10.^(sinclobe/10))/sum(10.^(omnilobe/10))
%   Comparison between sinc and omni.
% SEE ALSO
%  sinclobe

% by Sofia Mosesson 000626

if exist('gain', 'var')
 lobevector = ones(360,1)*gain;
else
 lobevector = ones(360,1)*3.0244;
end



