function res = beepsis_ndbin(data, dt, nbin, Gamma, methods)
% function res = beepsis_ndbin(data, dt, nbin,  Gamma, methods)
% 
% evaluates acceleration profile for the particle with inertia from its position
% (N degrees of freedom) using BEEPSIS approach
% force is assumed to be constant over one bin
% number of bins is input but could be decreased
%   maximal number of bins is defined in variable MAXBINS
% number of dagonals in the inversion of noise covariance matrix is given
%    by constant CINVTERMS
% 
% INPUT
%   data    matrix of positions (rows are time, columns degrees of
%              freedom), in meters
%   dt      time step between measured positions
%   nbin    number of bins
%             - empty     bin count automatically calculated
%             - 1 number  number of bins, same for all degrees of freedom
%             - vector of bin counts for, each number for 1 degree of
%                freedom
%   Gamma    damping coefficient           
%   methods  character string defing method of force  estimation. 
%            Can be any combination of following
%             'v'  inference for velocities (numerically calculated
%                  using central difference), 
%                  accorning to Eq. (E9)
%                  velocity calculated using backward difference
%             'p'  inference for positions and init velocity 
%                  velocity calculated using central difference
%                  according to Eq (E10)
%             'c'  inference for time correalted misfits in positions,
%                  according to Eqs (E4-E8)
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


MAXBIN = 100;
CINVTERMS = 25;


DIMS = size(data,2);
stddata = std(data, 'omitnan');
mindata = min(data);
maxdata = max(data);


%% bins count
if isempty(nbin)
    DX = (3.5*stddata)*size(data,1).^(-1/3);
    nbin = min(MAXBIN*ones(DIMS,1),floor((maxdata-mindata)/DX));
end
if numel(nbin) == 1
    nbin = nbin*ones(DIMS,1);
end

%% bin edges
xdata = cell(DIMS,1);
nbin = floor(nbin);
for idi = 1:DIMS
    xdata{idi} = linspace(mindata(idi), maxdata(idi), nbin(idi)+1);
    xdata{idi}(1)   = xdata{idi}(1)- eps(100);
    xdata{idi}(end) = xdata{idi}(end) + eps(100);
end

%% point indeces in bins, initial acc
GT    = Gamma*dt;
EGT   = exp(-GT);
M1EGT = -expm1(-GT);
GFP   = M1EGT/2/GT;

BINCOUNT  = prod(nbin);
DATACOUNT = size(data,1);

bincntr = deal(zeros(BINCOUNT,DIMS));
bincnt  = zeros(BINCOUNT,1);
binindex = zeros(size(data,1),1);
[forcev, forcep,  vdrift, a0]  = deal(zeros(BINCOUNT, DIMS));
[sumv, sump] = deal(0);
fullbins = 0;
fullbinsi = false(BINCOUNT,1);
for kk = 1:BINCOUNT
    ii = true(DATACOUNT,1);
    qi = cell(DIMS,1);
    [qi{:}] = ind2sub(nbin, kk);
    for idi = 1:DIMS 
        qq = qi{idi};
        ii = ii & (data(:,idi)>= xdata{idi}(qq)) & (data(:,idi) <=  xdata{idi}(qq+1) );
    end
    ii([1 end])  = false;
    
    if sum(ii) > 0
        fullbins = fullbins + 1;
        fullbinsi(kk) = true;
        binindex(ii) = fullbins;
        
        x0 = data(ii,:);
        xm = data(circshift(ii, -1),:);
        xp = data(circshift(ii, +1),:);
        
        vdrift(kk,:) = mean((xp-xm)/2/dt);
        a0(kk,:)     = mean((xp+xm-2*x0)/dt^2);
        
        
        xb = (xp-x0*(1+EGT) + xm*EGT)*Gamma/dt/M1EGT;
        forcev(kk,:)  = mean(xb,1);
        sumv         = sumv + var(xb, [],1);
        
        xb = (xp*(1-GFP)-x0 + xm*GFP)*Gamma^2/(GT-M1EGT);
        forcep(kk,:)  = mean(xb,1);
        sump         = sump + var(xb, [],1);
        
        bincntr(kk,:) = mean(x0,1);
        bincnt(kk)    = bincnt(kk) + sum(ii);
    else
        forcev(kk,:)  = NaN;
        vdrift(kk,:) = NaN;
        a0(kk,:)     = NaN;
        
        forcep(kk)  = NaN;
        
        for idi = 1:DIMS
            qq = qi{idi};
            bincntr(kk,idi) = mean(xdata{idi}(qq:qq+1));
        end
        
    end
    
end


%% data init
x0  = data(2:end-1,:);
xp1 = data(3:end,:);
xm1 = data(1:end-2,:);
binindex = binindex(2:end-1);

%% process
nbin = nbin(:)';

res.bincenter = reshape(bincntr, [nbin, DIMS]);
res.binedges  = xdata;
res.bincount  = reshape(bincnt, nbin);
res.vdrift    = reshape(vdrift, [nbin, DIMS]);  
res.a0        = reshape(a0, [nbin, DIMS]);  


if any(methods == 'v') % Eq (B1)
    res.forcev = reshape(forcev,  [nbin, DIMS]);
end

if any(methods == 'p') % Eq (B2)
    res.forcep = reshape(forcep,  [nbin, DIMS]);
end
 
if any(methods == 'c') % Eqs (7-11, 19)
    p = zeros(fullbins,DIMS);
    for idi = 1:DIMS
        pl = calc_force(xp1(:,idi), x0(:,idi), xm1(:,idi), binindex, Gamma, dt, CINVTERMS);
        p(:, idi) = pl;
    end
    
    pp = NaN(prod(nbin), DIMS);
    pp(fullbinsi,:) = p;
    pp = reshape(pp, [nbin, DIMS]);
    res.force = pp;
end


end


function [par, err, CM, fval] = calc_force(xp, x0, xm, binindex, GM, dt, CINVTERMS)
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
    
    
    for kcol = 1:CINVTERMS
        bj = binindex(kcol);
        mc = xb(kcol);
        for krow = 1:kcol+CINVTERMS
            bi = binindex(krow);
            Dij = (-RP).^(-abs(krow-kcol))/RD * (1-RR.^(min(krow, kcol)));
            PM(bi) = PM(bi) + Dij*mc;
            IM(bi, bj) = IM(bi, bj) + Dij;
        end
    end
    
    for kcol = N-CINVTERMS+1:N
        bj = binindex(kcol);
        mc = xb(kcol);
        for krow = kcol-CINVTERMS:N
            bi = binindex(krow);
            Dij = (-RP).^(-abs(krow-kcol))/RD * (1-RR.^(N+1-max(krow, kcol)));
            PM(bi) = PM(bi) + Dij*mc;
            IM(bi, bj) = IM(bi, bj) + Dij;
        end
    end
    
    
    kk  = -CINVTERMS:CINVTERMS;
    Dij = (-1).^(-abs(kk)).*RP.^(-abs(kk))/RD;
    for kcol = CINVTERMS+1:N-CINVTERMS
        bj = binindex(kcol);
        mc = xb(kcol);
        id  = 1;
        for krow = kcol-CINVTERMS:kcol+CINVTERMS
            bi = binindex(krow);
            PM(bi) = PM(bi) + Dij(id)*mc;
            IM(bi, bj) = IM(bi, bj) + Dij(id);
            id = id + 1;
        end
    end
    
    
    par  = IM\PM;
    err  = NaN(size(par));
    CM   = NaN(numel(par));
    fval = 0;
end

