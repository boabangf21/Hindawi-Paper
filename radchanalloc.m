function [k, obkn] = radchanalloc(b, obk, bdes)
% DESCRIPTION [k ,obkn] = radchanalloc(b, obk, bdes) 
%  Radio Channel Allocation. 
% INPUT
%  b --      Bases where channels are requested.
%  obk --    Binary matrix, base by channel.
% OUTPUT
%  k --      Channel allocation vector, 
%            related to b and with same size.
%            NaNs where allocation failed.
%  obkn --   The new obk matrix when allocation is made.
%  bdes --   Desired base station.
% TRY 
%  [k, obkn] = radchanalloc([1], [1 1; 1 1])
%  [k, obkn] = radchanalloc([2], [1 1; 1 1])
%  [k, obkn] = radchanalloc([2 2 2], [1 1; 1 1])
%  [k, obkn] = radchanalloc(nan, [1 1; 1 1])
%  radchanalloc(1,   1    ) => [1,   0     ]
%  radchanalloc(1,   0    ) => [nan, 0     ]
%  radchanalloc(nan, 1    ) => [nan, 1     ]
%  radchanalloc(nan, 0    ) => [nan, 0     ]

% by Magnus Almgren 990401

obkn = obk; 
bv = b(~isnan(b));
k  = nans(size(b));
kv = nans(size(bv));
if length(bv(:))>0&length(obk(:))>0
 req = full(sparse((1:size(bv(:),1))',bv,1));
 
 [aa, bb] = padsiz(0, req, obk'); % Make to same size. Pad with 0.
 ra = sqrt(eps)*rand(size(aa));
 rb = sqrt(eps)*rand(size(aa));
 [sa,ia] = sort(ra-aa,1); % Sort along first dim.
 [sb,ib] = sort(rb-bb,1);
 s = (-sa>0.5)&(-sb>0.5); % Look for common space.
 
 bind = index(ib,1:size(bb,2),size(bb));
 bb(bind(s)) = 0;
 obkn = bb';
 obkn = obkn(1:size(obk,1),1:size(obk,2)) > 0.5;
 
 kv(ia(s)) = ib(s);
end
k(~isnan(b)) = kv;















