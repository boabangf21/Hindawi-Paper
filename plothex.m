function plothex(xyb, fib, chgroup)
% DESCRIPTION plothex(xyb, fib, chgroup)
%  Plots the homogeneous hexagon cell pattern with BS identities. If chgroup
%  is specified the hexagons are filled with a color.
% INPUT
%  xyb --      Base station coordinates, given as a complex row vector.
%  fib --      The base station antenna vector, given as a complex row
%              vector pointing from base to middle of cell.
%  chgroup --  Row vector. Contains the channel group number of the cell.
%              This number is used to set the color of the cell.
% TRY
%  plothex(0,0)
%  plothex([0 1],[0 0])
%  plothex(irandn(3,1),irandn(3,1),(1:3))
% SEE ALSO
%  crecells

% By Fredric Kronestedt 960212 modified by MA 000522 OQ 000523

% Preserve the hold state
holdState=ishold;

% Calculate the cell radius
if all(fib==0)       
 % omni cell (fib==0)
 b2b = sort(abs(xyb-xyb(1))); % shortest dist between cellcenters
 if length(b2b)>1
  cellradius = mplus(fib,b2b(2)*i/sqrt(3));
 else
  cellradius = 1000i; % the case with one omnibase only
 end
 tpos = cellradius(1)/4; % Place text slightly off the cell center
else
 % Directional antennas
 cellradius = fib;
 tpos = 0; % Place text in the middle of the cell.
end

% Plot the hexagonal pattern
% Each column of 'points' contain the corner points of each cell
points=mplus(xyb,fib,exp(i*pi*(0:6)/3).'*cellradius);
if(exist('chgroup','var'))
 % Channel groups are specified. Draw filled hexagons
 cmap=rgb2hsv(hsv(max(chgroup)));
 cmap=mprod(cmap,[1 0.3 1]); % Create unsaturated colors and set as colormap
 colormap(hsv2rgb(cmap));
 % We must have a for loop to avoid problems when we have 3 cells
 for j=1:size(chgroup,2) 
  fill(real(points(:,j)),imag(points(:,j)),chgroup(j));
  hold on; % Hold the image so that the text and hexagons show up at the same plot.
 end
else
 % No colors are specified. Just plot the borders
 plot(real(points),imag(points),'k');
 hold on; % Hold the image so that the text and hexagons show up at the same plot.
end


% Plot base positions
plot(real(xyb),imag(xyb),'ro') 
axis('equal')

% Plot BS id
ts=get(gcf,'DefaultTextFontSize'); % Preserve font size
set(gcf,'DefaultTextFontSize',7);  % Set font size to 7
pos = xyb + fib + tpos;            % BS text pos
for s=1:length(xyb)
 text(real(pos(s)),imag(pos(s)),num2str(s),'HorizontalAlignment','center')
end
set(gcf,'DefaultTextFontSize',ts); % Restore old font size

  
% Restore the hold state (if it was off, it has been set to on earlier)
if(~holdState)
 hold off;
end