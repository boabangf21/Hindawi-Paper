function lobevector = sinclobe
% DESCRIPTION lobevector = sinclobe
%  Creates an antenna gain lobe vector. The antenna lobe is approximated
%  with a log(sinc) within the sector and with -30 dB outside the sector.
% OUTPUT
%  lobevector --  A row vector containing the gain in each 1 degree direction.
%                 The length of the vector is 360 where last index corresponds 
%                 to zero degrees and first index to 1 degree.
% TRY
%  plot(sinclobe), xlabel('[degrees]'), ylabel('gain [dB]')
% SEE ALSO
%  pathgain

% by Magnus Almgren 000512

lw=2.256;   % Fitted to give -10.53 db gain @ 60 degrees <=> 35*log10(0.5).
fi = (1:75)/360*2*pi; 
lobevector = max(-30,10/log(10)*log((sin(fi*lw)./(fi*lw)).^2));
lobevector = [lobevector zeros(1,104)-30];
lobevector = [lobevector -30 fliplr(lobevector) 0]';
lobevector = lobevector+10; % 10 dB gain in main direction.




