function par = setpar
% DESCRIPTION par = setpar
%  Sets the structured parameter par that is the input to runef.
%  This file is intended to use as a template. Rename and change
%  its content to your desire.
% OUTPUT      
%  par --  input parameter to rune
% TRY
%  setpar
% SEE ALSO
%  runef

% by Magnus Almgren 000530

 par.cellradius = 1000; % cell radius [m]
 par.km = 1;            % km^2+lm^2+km*lm => the number of sites
 par.lm = 1;            % related to km above
 par.sps = 3;           % Number of sectors per site
 par.ncluster = 3;      % number of clusters equals ncluster^2
 par.kpc = 2;           % number of allocated channels per cell
 par.gainconst = 21;    % gain at 1 meter distance 
 par.noise = -118;      % thermal noise floor -118 dBm
 par.alpha = 3.5;       % distance attenuation coefficient
 par.sigma = 6;         % standard deviation for the log-normal fading in dB
 par.raa = 0.5;         % lognormal correlation down link (typical 0.5)
 par.corrdist = 110;    % lognormal fading correlation distance [m]
 par.offtraf = 2;       % average number of offered calls to a cell     
 par.mht = 20;          % mean holding time [seconds]  
 par.amean = 1;         % mean acceleration [m/s/s]
 par.vmean = 14;        % mean speed  [m/s]
 par.pinit = 20;        % init power set to each new link
 par.sirmin = 10;       % C/I level under which a call is dropped  
 par.homargin = 3;      % gain margin between two bases used at Hand Off 
 par.dt = 1;            % time interval in the simulation loop
 par.maxiter = 10;      % number of iterations in the main loop 
 par.seed = 1;          % seed to all random sequencies in the simulation




