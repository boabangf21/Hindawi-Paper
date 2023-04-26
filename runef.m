function [res, par, sys] = runef(par,sta,sys) 
% DESCRIPTION [res, par, sys] = runef(par,sta,sys)
%   The basic rune dynamic function.
% INPUT  (all inputs are optional)
%  par --            BASIC SIMULATION PARAMETERS
%  par.cellradius -- cell radius [m]
%  par.km --         km^2+lm^2+km*lm => the number of sites
%  par.lm --         related to km above
%  par.sps --        number of sectors per site
%  par.ncluster --   number of clusters equals ncluster^2
%  par.kpc --        number of allocated channels per cell
%  par.gainconst --  gain at 1 meter distance [dB] 
%  par.noise --      thermal noise floor [dBm]
%  par.alpha --      distance attenuation coefficient
%  par.sigma --      standard deviation for the lognormal fading [dB]
%  par.raa --        lognormal correlation down link (typical 0.5)
%  par.corrdist --   lognormal fading correlation distance [m]
%  par.offtraf --    average number of offered calls to a cell [Erlang/cell]     
%  par.mht --        mean holding time [seconds]  
%  par.amean --      acceleration [m/s/s]
%  par.vmean --      speed  [m/s]
%  par.pinit --      init power set to each new link [dBm]
%  par.sirmin --     C/I level under which a call is dropped [dB] 
%  par.homargin --   gain margin between two bases used at Hand Off [dB]
%  par.dt --         time interval in the simulation loop [s]
%  par.maxiter --    number of iterations in the main loop 
%  par.seed --       seed to all random sequencies in the simulation
%  sta --            STATE
%  sta.xym --        position in complex form [m]
%  sta.xyv --        speed imag <=> north real <=> east [m/s] 
%  sta.m --          mobile identity number
%  sta.b --          base station index
%  sta.k --          channel number
%  sta.pul --        transmitted power up link  [dBm]
%  sta.cul --        received carrier power up link [dBm]
%  sta.iul --        interference power up link  [dBm]
%  sta.sirul --      signal to interference ratio up link [dB]
%  sta.pdl --        transmitted power down link  [dBm]
%  sta.cdl --        carrier downlink [dBm]
%  sta.idl --        interference down link [dBm]
%  sta.sirdl --      signal to interference down link [dB]
%  sta.mtop --       max mobile number used sofar
%  sta.obk --        the allocation of channels
%  sta.seed --       the seed before next iteration
%  sys --            ALL INTERMEDIATE VALUES NEEDED TO SIMULATE
%  sys.xyb --        base positions [m]
%  sys.fib --        cell center vector [m]
%  sys.rombvec --    stem vectors [m]
%  sys.lobevector -- antenna gain diagram [dB]
%  sys.iniobk --     channel plan
%  sys.lognmap --    lognormal map [dB]
%  sys.mapvec --     lognormal map vectors [m]
% OUTPUT
%  res --          cellarray of sta collected for each iteration of the simulation
%  par --          same as input otherwise created within the function
%  sys --          same as input otherwise created within the function  
% TRY 
%  [res, par, sys] = runef 

% by Magnus Almgren 000517 	

% set simulation parameter par if not present as an input
if ~exist('par','var')
 par = setpar; % default parameter setting 
end

% Create the sys variable if not present as an input.
if ~exist('sys','var')
 % generate base station position and directions 
 [sys.xyb, sys.fib, sys.rombvec] = crecells(par.cellradius,par.sps,par.km,par.lm,par.ncluster);
 % Antenna gain for all directions, size == [360 1].
 if all(abs(sys.fib)>0)
  sys.lobevector = sinclobe;    
 else 
  sys.lobevector = omnilobe;
 end
 
 % Create a channel plan for the system.
 % Number of channels that is used in each cluster.   
 nk = length(sys.xyb)/par.ncluster.^2*par.kpc;    	
 sys.iniobk   = crechanplan(length(sys.xyb),nk,par.ncluster);  % Allocate channels to cells
 
 % Create a lognormal map.
 if par.sigma > 0  % Is a lognormal map needed (takes a few seconds to generate).
  % The lognormal map is dependent on the seed.
  oseed = setseed(par.seed);  % Set seed of pseudo random generator for the map.
  [sys.lognmap, sys.mapvec] = crelognmap(sys.xyb, sys.rombvec, par.corrdist); 
  setseed(oseed);  % Restore seed to original value.
 else
  sys.lognmap = 0; % Give fake arguments to pathgain,
  sys.mapvec = 0;  % the values doesn't matter anyway.
 end
end


% init of state variable
% The variables below are altered in the for loop and save after
% each iteration
if ~exist('sta','var')
 e = zeros(0,1);
 sta.xym   = e; 
 sta.xyv   = e; 
 sta.m     = e; 
 sta.b     = e; 
 sta.k     = e;
 sta.pul   = e; 
 sta.cul   = e; 
 sta.iul   = e; 
 sta.sirul = e;
 sta.pdl   = e; 
 sta.cdl   = e; 
 sta.idl   = e; 
 sta.sirdl = e;
 sta.mtop  = 0; 
 sta.obk   = sys.iniobk;    
 sta.seed  = par.seed;
end

oseed=setseed(sta.seed); % Set seed in random generators. SM 000713

% The simulation loop.
for iter = 1:par.maxiter
 
 % Terminate some calls and drop calls with bad quality. 
 terones = ~isnan(sta.k)&(rand(size(sta.xym)) < 1-exp(-mdiv(par.dt,par.mht))); % Natural terminated calls.
 drones = (min([sta.sirul sta.sirdl],[],2) < par.sirmin);  % quality dropping
 
 % Terminate by setting b and k to NaN and obk to 1
 [sta.b,sta.k,sta.obk] = terminate(sta.b,sta.k,sta.obk,drones|terones);
  
 % Make a realisation of new users.
 nmob = mrequest(par.offtraf * length(sys.xyb), mdiv(par.dt,par.mht));
 nt = nans(nmob,1);  % NaN vector used to concatenate with below.
 mn = sta.mtop+(1:nmob)'; % new mobile id numbers
 sta.mtop = sta.mtop+nmob; % highest mobile id so far 
 
 % Clean variables. 
 e = isfinite(sta.k); % Keep only live calls without channel
 
 % Refresh vector structure by removing released calls and adding new ones.
 sta.xym=[sta.xym(e); nt]; 
 sta.xyv=[sta.xyv(e); nt]; 
 sta.m=[sta.m(e); mn]; 
 sta.b=[sta.b(e); nt]; 
 sta.k=[sta.k(e); nt]; 
 sta.pul=[sta.pul(e); nt]; 
 sta.cul=[sta.cul(e); nt]; 
 sta.iul=[sta.iul(e); nt]; 
 sta.sirul=[sta.sirul(e); nt];
 sta.pdl=[sta.pdl(e); nt]; 
 sta.cdl=[sta.cdl(e); nt]; 
 sta.idl=[sta.idl(e); nt]; 
 sta.sirdl=[sta.sirdl(e); nt];
 
 % Move old users and initiate a position to newcomers.
 [sta.xym,sta.xyv] = mobmove(sta.xym,sta.xyv,par.amean,par.vmean,par.dt,sys.rombvec);

 % Calculate the gain matrix.
 sta.g = pathgain(sta.xym,  sys.xyb, sys.fib, sys.rombvec, ...
  par.gainconst, par.alpha, par.sigma, par.raa, ...
  sys.lobevector, sys.lognmap, sys.mapvec);

 % Perform handoff & Allocate new mobiles to a channel if possible.
 [sta.b,sta.k,sta.obk] = handoff(sta.b,sta.k,sta.g,sta.obk,par.homargin);
 
 % Try to assign channels to new calls.
 [sta.b,sta.k,sta.obk] = assign(sta.b,sta.k,sta.g,sta.obk,par.homargin);
 
 % Set power to new users.
 sta.pul(isnan(sta.pul)) = par.pinit; 
 sta.pdl(isnan(sta.pdl)) = par.pinit;
 
 % Reset power of terminated users.
 sta.pul(isnan(sta.k)) = nan; 
 sta.pdl(isnan(sta.k)) = nan;
 
 % Perform the transmission from the sending to the receiving side,
 % all in dB or dBm.
 [sta.cul, sta.iul, sta.sirul] = transmitul(sta.b, sta.k, sta.pul, sta.g, par.noise); % up link
 [sta.cdl, sta.idl, sta.sirdl] = transmitdl(sta.b, sta.k, sta.pdl, sta.g, par.noise); % down link
 
 % Collect results.
 sta.seed = setseed;  % the seed before next iteration
 res{iter} = sta;     % Save current state.
end

setseed(oseed); % Restore seed to original value. Sofia Mosesson 000713