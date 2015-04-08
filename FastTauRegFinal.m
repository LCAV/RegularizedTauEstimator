function result = FastTauRegFinal(x, y,lambda, control)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes regularized tau estimate of regression
%
% tau-estimate is tuned to have 95% efficiency, and 50% bdp,
% using Optimal rho-function

%
% INPUT:
% x : mixing matrix
% y : vector of measurements
% control.N : number of random initial solutions (e.g. 500)
%        .k : number of initial IRWLS steps on each candidate (e.g. 2)
%        .t : number of best candidates to fully improve (e.g. 5)
%        .r : number of iterations in scale approximation in case approx=1 (e.g. 2)
%        .approx : if 0, fully compute S-scale, otherwise approximate
%                   (approx=0 is recommended)
%       
%
% OUTPUT:
% result.beta : the tau-estimate of regression
%       .scale : the tau-estimate of scale
%
% Author : Marta Martinez-Camara,LCAV-EPFL, 2015
%          marta.martinez-camara@epfl.ch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rand('state', 1); % this would set the seed for random subsampling
%                           (and make result reproducable; but be careful in simulations!)

if control.t < 1
    error('parameter t should be at least 1')
end

N = control.N;
k = control.k;
bestr = control.t;
scalesteps = control.r;
approx = control.approx;

RWLStol = 1e-11; % tolerance for IRWLS convergence
Mscaletol = 1e-5; % tolerance for M-scale iteration convergence

[n,p] = size(x);
c1 = 0.4046; % c parameter for rho_1 function
b = .5;      % b parameter for M-scale estimation
c2 = 1.0900; % c parameter for rho_2 function


bestbetas = zeros(bestr,p); % storing the t best candidates
bestscales = 1e20 * ones(bestr,1); % storing the scale of the t best candidates
besttauscales = 1e20 * ones(bestr,1); % storing tau-scale of the t best candidates
worstind = 1; % the index of the worst one in these best ones
worsttau = 1e20; % the worst tau of the bestr best ones
worsts = 1e20; % the S-scale corresponding to the worst tau
worstres = y;

for i=1:N
    % draw a non-singular random subsample
    singular = 1; itertest=1;
    while (singular == 1 && itertest < 100)
        indices = randperm(n);
        xs = x(indices(1:p),:);
        ys = y(indices(1:p));
        if rank(xs) == p
            bbeta = xs\ys;% creating a random initial solution
            singular = 0;
        else
            itertest = itertest+1;
        end
    end
    if itertest==100
        error('too many degenerate subsamples')
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % do k initial IRLS-steps
    if k>0
        tmp = IRLSiteration(x, y, bbeta, 0, k, RWLStol, b, c1, c2,lambda);
        betarw = tmp.betarw;
        resrw = y - x * betarw;
        scalerw = tmp.scalerw;
    else
        betarw = bbeta;
        resrw = y - x * betarw;
        scalerw = median(abs(resrw))/.6745;
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check whether new subsample is one of best t 

    if ~approx % compute actual scale, but use tau-conditions!
        
            news = solveMscale(resrw, scalerw, Mscaletol, b, c1);
            newtau = news * sqrt(mean(rhoOptfun(resrw/news,c2)))+lambda*(betarw'*betarw);
            if newtau < worsttau
                besttauscales(worstind) = newtau;
                bestscales(worstind) = news;
                bestbetas(worstind,:) = betarw';
                [worsttau, worstind] = max(besttauscales);
                worsts = bestscales(worstind);
                worstres = y - x * bestbetas(worstind,:)';
            end
        
    
    else % or just compute approximations (and don't bother with the conditions)
        news = scalerw;
        for kstep = 1:scalesteps
            news = sqrt( news^2 * mean( rhoOptfun(resrw/news,c1) ) / b );
        end
        newtau = news * sqrt(mean(rhoOptfun(resrw/news,c2)))+lambda*(betarw'*betarw);
        if newtau < worsttau
            besttauscales(worstind) = newtau;
            bestscales(worstind) = news;
            bestbetas(worstind,:) = betarw';
            [worsttau, worstind] = max(besttauscales);
            worsts = bestscales(worstind);
            worstres = y - x * bestbetas(worstind,:)';
        end
    end
end

%%%%%%%%%%%%%%%%
% First k IRLS iterations for this initial solution

if k>0
    tmp = IRLSiteration(x, y, bbeta, 0, k, RWLStol, b, c1, c2,lambda);
    betarw = tmp.betarw;
    resrw = y - x * betarw;
    scalerw = tmp.scalerw;
else
    betarw = bbeta;
    resrw = y - x * betarw;
    scalerw = median(abs(resrw))/.6745;
end
if ~approx
    
        news = solveMscale(resrw, scalerw, Mscaletol, b, c1);
        newtau = news * sqrt(mean(rhoOptfun(resrw/news,c2))) + lambda*(betarw'*betarw);
        if newtau < worsttau
            bestscales(worstind) = news;
            bestbetas(worstind,:) = betarw';
        end
    
else
    news = scalerw;
    for kstep = 1:scalesteps
        news = sqrt( news^2 * mean( rhoOptfun(resrw/news,c1) ) / b );
    end
    newtau = news * sqrt(mean(rhoOptfun(resrw/news,c2))) + lambda*(betarw'*betarw);
    if newtau < worsttau
        bestscales(worstind) = news;
        bestbetas(worstind,:) = betarw';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterate the t best candidates until convergence

superbesttauscale = 1e20;
for i = bestr:-1:1
    tmp = IRLSiteration(x, y, bestbetas(i,:)', bestscales(i), 1000, RWLStol, b, c1, c2,lambda);
    resrw = y - x * tmp.betarw;
    tauscalerw = tmp.scalerw * sqrt(mean(rhoOptfun(resrw/tmp.scalerw,c2))) + lambda*(tmp.betarw'*tmp.betarw);
    if tauscalerw < superbesttauscale
        superbesttauscale = tauscalerw;
        superbestbeta = tmp.betarw;
        superbestscale = tmp.scalerw;
    end
end

superbestscale = solveMscale(y-x*superbestbeta, superbestscale, Mscaletol, b, c1);
superbesttauscale = superbestscale * sqrt(mean(rhoOptfun((y-x*superbestbeta)/superbestscale,c2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add TIK inital as an initial solution

betaLS = (x'*x + (lambda)*eye(p))\x'*y; % Tikhonov solution
resLS = y - x * betaLS;
scaleLS = median(abs(resLS))/.6745;
news = solveMscale(resLS, scaleLS, Mscaletol, b, c1);
newtau = news * sqrt(mean(rhoOptfun(resLS/news,c2)))+lambda*(betaLS'*betaLS);
if newtau < superbesttauscale
     %superbestscale = news;
     superbestbeta = betaLS;
     superbesttauscale = newtau;
end   


result.beta = superbestbeta;
result.scale = superbesttauscale 

%--------------------------------------------------------------------------

function result = IRLSiteration(x, y, initialbeta, initialscale, maxiter, tol, b, c1, c2,lambda)

% approximate IRLS iteration; pass maxiter=500, say, if convergence is desired
% e.g. tol = 1e-11

[n,p]=size(x);
res = y - x * initialbeta;
if (initialscale == 0)
    scale = median(abs(res))/.6745;
else
    scale = initialscale;
end
beta = initialbeta;

betadiff=2*tol;
iter = 0;
while (betadiff > tol) && (iter < maxiter)
    scale = sqrt( scale^2 * mean( rhoOptfun(res/scale,c1) ) / b );
    scaledres = res/scale;
    Wn_numer = sum(WnumerOptfun(scaledres,c2));
    Wn_denom = sum(psiOptxfun(scaledres,c1));
    Wn = Wn_numer/Wn_denom;
    weights = (Wn*fw(scaledres,c1)+fw(scaledres,c2))/2*n; 
    sqweights = weights.^(1/2);
    sqW = sqweights * ones(1,p);
    xw = x .* sqW;
    yw = y .* sqweights;
    newbeta=(xw'*xw + (lambda)*eye(p))\xw'*yw;
    if (any(isnan(newbeta)))
        newbeta = initialbeta;
        scale = initialscale;
        break
    end

    betadiff = norm(beta - newbeta)/sqrt(p);
    res = y - x * newbeta;
    beta = newbeta;
    iter = iter + 1;
end

result.betarw = newbeta;
result.scalerw = scale;
result.iters = iter;

%--------------------------------------------------------------------------

function sc = solveMscale(x, initialscale, tol, b, c)

% M-estimator of scale using the Optimal Rho function.

maxiter = 100;
if (initialscale==0)
    s=median(abs(x))/.6745;
else
    s=initialscale;
end

rhoold = mean(rhoOptfun(x/s,c)) - b ;
iter = 0;
while (abs(rhoold) > tol) && (iter < maxiter)
    delta = rhoold / mean( psiOptxfun(x/s,c) ) / s;
    isqu = 1; ok = 0;
    while  (isqu < 30 && ok~=1)
        rhonew = mean( rhoOptfun(x/(s+delta),c)) - b;
        if abs(rhonew) < abs(rhoold)
            s = s + delta; ok = 1;
        else
            delta = delta/2 ; isqu = isqu + 1 ;
        end
    end
    if isqu==30
        maxiter = iter; % we tell it to stop, but we keep the iter for info
    end
    rhoold = rhonew;
    iter = iter + 1;
end
sc = abs(s);

%--------------------------------------------------------------------------

function psi = fw(x,c)
% weight function = psi(x)/x

inds = 1:length(x);
tmp = zeros(length(x),1);
tmp(abs(x) <= 2*c) = 1 / (3.25*c^2);
fullinds = inds((abs(x) > 2*c)&(abs(x) <= 3*c));
fullx = x(fullinds);
tmp(fullinds) = (-1.944  / c^2 + 1.728 * fullx.^2 / c^4 - 0.312 * fullx.^4 / c^6 + 0.016 * fullx.^6 / c^8) / 3.25;
psi = tmp;


%--------------------------------------------------------------------------

function rho=WnumerOptfun(x, c)

% Computes function in numerator of W

inds = 1:length(x);
tmp = zeros(length(x),1);
fullinds = inds((abs(x) > 2*c)&(abs(x) <= 3*c));
fullx = x(fullinds);
tmp(fullinds) = (3.584 - 0.864 * fullx.^4 / c^4 + 0.208 * fullx.^6 / c^6 - 0.012 * fullx.^8 / c^8) / 3.25;
tmp(abs(x) > 3*c) = 2;
rho = tmp;



%--------------------------------------------------------------------------

function rho=rhoOptfun(x, c)

% Computes Optimal rho function

inds = 1:length(x);
tmp = ones(length(x),1);
inds1 = inds(abs(x) <= 2*c);
tmp(inds1) = x(inds1).^2 / 2 / (3.25*c^2);
fullinds = inds((abs(x) > 2*c)&(abs(x) <= 3*c));
fullx = x(fullinds);
tmp(fullinds) = (1.792 - 0.972 * fullx.^2 / c^2 + 0.432 * fullx.^4 / c^4 - 0.052 * fullx.^6 / c^6 + 0.002 * fullx.^8 / c^8) / 3.25;
rho = tmp;


% --------------------------------------------------------------------

function psi=psiOptxfun(x,c)

% Computes Optimal Rho function's first derivative times x

inds = 1:length(x);
tmp = zeros(length(x),1);
inds1 = inds(abs(x) <= 2*c);
tmp(inds1) = x(inds1).^2 / (3.25*c^2);
fullinds = inds((abs(x) > 2*c)&(abs(x) <= 3*c));
fullx = x(fullinds);
tmp(fullinds) = (-1.944 * fullx.^2 / c^2 + 1.728 * fullx.^4 / c^4 - 0.312 * fullx.^6 / c^6 + 0.016 * fullx.^8 / c^8) / 3.25;
psi = tmp;



