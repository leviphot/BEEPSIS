function res = beepsis_bin(data, dt, bins,  Gamma, methods)
% function res = beepsis_bin(data, dt, xdata,  Gamma, methods)
% 
% evaluates acceleration profile for the particle with inertia from its position
% using BEEPSIS approach
% force is assumed to be constant over one bin
% number of bins is input but could be decreased
%   maximal number of bins is defined in variable MAXBINS
%   minimal number of points in 1 bin is NMIN
%     if number of positions is lower, adjacent bins are merged
% number of dagonals in the inversion of noise covariance matrix is given
%    by constant CINVTERMS
% 
% INPUT
%   data     vector of positions (1D), should be in meters
%   dt       time step between measured positions
%   bins     number of bins/bin edges
%             - empty     bin count automatically calculated
%             - 1 number  maximal number of bins
%             - vector of bin edges  
%   Gamma    damping coefficient           
%   methods  character string defing method of force  estimation. 
%            Can be any combination of following
%             'v'  inference for velocities (numerically calculated
%                  using central difference), 
%                  accorning to Eq. (E9)
%                  velocity calculated using backward difference
%             'p'  inference for positions and init velocity based on
%                  Chandrasekar work
%                  velocity calculated using central difference
%                  according to Eq (E10)
%             'c'  inference for time correalted misfits in positions,
%                  according to Eqs (E7-E8)
%                  
% OUTPUT
%   res      structure with results
%            bincenter  center of a given bin
%            binegdes   edges of bins
%            forcev     force calculated using Eq. (E9)
%            forcep     force calculated using Eq. (E10)
%            force      force calcuated according Eqs (E4-E8)
%            vdrift     drift velocity averaged in each bin
%
% This file is part of the BEEPSIS toolbox.
% See LICENSE.md for information about using/distributing this file.

NMIN      = 20;
MAXBIN    = 100;
CINVTERMS = 25;

data = data(:);
stddata = std(data, 'omitnan');
mindata = min(data);
maxdata = max(data);



%% bins count
if isempty(bins)
    DX = (3.5*stddata)*numel(data).^(-1/3);
    bins = min(MAXBIN,floor((maxdata-mindata)/DX));
end

%% bin edges
if numel(bins) == 1
    bins = floor(bins);
    bins = linspace(mindata, maxdata, bins+1);
    bins(1)   = bins(1)   -  10*eps(1);
    bins(end) = bins(end) +  10*eps(1);
    nbin = zeros(numel(bins)-1,1);
    for kk = 1:numel(nbin)
        ii = data>bins(kk) & data<=bins(kk+1);
        nbin(kk) = sum(ii);
    end
    while true        
        [mb, imin] = min(nbin);
        if (mb >= NMIN) || numel(bins)<=10
            break
        end
        
        if imin == 1
            ijoin = 2;
        elseif imin == numel(nbin)
            ijoin = numel(nbin)-1;
        elseif nbin(imin - 1) < nbin(imin + 1)
            ijoin = imin - 1;
        else
            ijoin = imin + 1;
        end
        if ijoin < imin
            imin = ijoin;
        end
        bins = [bins(1:imin) bins(imin+2:end)];
        nbin = [nbin(1:imin-1); (nbin(imin)+nbin(imin+1)); nbin(imin+2:end)];        
    end 
        
end

%% point indeces in bins, forces using Eqs (B1 and B2)
GT    = Gamma*dt;
EGT   = exp(-GT);
M1EGT = -expm1(-GT);
GFP   = M1EGT/2/GT;


bincntr = deal(zeros(numel(bins)-1,1));
binindex = zeros(numel(data),1);
[forcev, forcep]  = deal(zeros(1,numel(bins)-1));
for kk = 1:numel(bins)-1
    ii = (data>= bins(kk)) & (data <  bins(kk+1) );
    ii([1 end])  = false;
    
    binindex(ii) = kk;
    
    x0 = (data(ii));
    xm = (data(circshift(ii, -1)));
    xp = (data(circshift(ii, +1)));
    
    xb = (xp-x0*(1+EGT) + xm*EGT)*Gamma/dt/M1EGT;
    forcev(kk)  = mean(xb);    
    
    xb = (xp*(1-GFP)-x0 + xm*GFP)*Gamma^2/(GT-M1EGT);    
    forcep(kk)  = mean(xb);
    
    bincntr(kk) = mean(x0);
    
end


%% data init
x0  = data(2:end-1);
xp1 = data(3:end);
xm1 = data(1:end-2);
binindex = binindex(2:end-1);

if (sum(binindex==0) ~= 0)
    error('BEEPSISBIN:point_out_of_bins', 'Some postion(s) was not assigned to a bin');
end

%% process


res.bincenter = bincntr;
res.binedges  = bins;
res.bincount  = nbin;


if any(methods == 'v') % Eq (B1)
    res.forcev = forcev;
end

if any(methods == 'p') % Eq (B2)
    res.forcep = forcep;
end

if any(methods == 'c') % Eqs (7-11, 19)
    res.force = calc_force(xp1, x0, xm1, binindex, Gamma, dt, CINVTERMS);
end


end


function force = calc_force(xp, x0, xm, binindex, GM, dt, NQ)
    GT    = GM*dt;
    EGT   = exp(-GT);
    M1EGT = -expm1(-GT);
    EGT2 = EGT^2;

    if GT < 1e-4
        D1  = 4/3 - 4/3*GT + 4/5*GT^2;
        D2  = 1/3 - 1/3*GT + 11/60*GT^2;
    else
        D1 = (2*GT-2  + (2*GT+2)*EGT2)/GT^3;
        D2 = (1-EGT2-2*GT*EGT)/GT^3;
    end
    xb = (xp-x0*(1+EGT) + xm*EGT)*GM/dt/M1EGT;

    xcheb = D1/2/D2;
    RP    = xcheb + sqrt(xcheb^2-1);
    RM    = xcheb - sqrt(xcheb^2-1);
    RR    = RM/RP;
    RD    = RP-RM;

    BC = max(binindex);
    N  = numel(x0);
    
    PM = zeros(BC,1);
    IM = zeros(BC);
    
    
    for kcol = 1:NQ
        bj = binindex(kcol);
        mc = xb(kcol);
        for krow = 1:kcol+NQ
            bi = binindex(krow);
            Dij = (-RP).^(-abs(krow-kcol))/RD * (1-RR.^(min(krow, kcol)));
            PM(bi) = PM(bi) + Dij*mc;
            IM(bi, bj) = IM(bi, bj) + Dij;
        end
    end
    
    for kcol = N-NQ+1:N
        bj = binindex(kcol);
        mc = xb(kcol);
        for krow = kcol-NQ:N
            bi = binindex(krow);
            Dij = (-RP).^(-abs(krow-kcol))/RD * (1-RR.^(N+1-max(krow, kcol)));
            PM(bi) = PM(bi) + Dij*mc;
            IM(bi, bj) = IM(bi, bj) + Dij;
        end
    end
    
    
    kk  = -NQ:NQ;
    Dij = (-1).^(-abs(kk)).*RP.^(-abs(kk))/RD;
    for kcol = NQ+1:N-NQ
        bj = binindex(kcol);
        mc = xb(kcol);
        id  = 1;
        for krow = kcol-NQ:kcol+NQ
            bi = binindex(krow);
            PM(bi) = PM(bi) + Dij(id)*mc;
            IM(bi, bj) = IM(bi, bj) + Dij(id);
            id = id + 1;
        end
    end
    
    
    force  = IM\PM;
end


