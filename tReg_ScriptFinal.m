%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to generate Figure 2 from the paper 'A new robust 
% and efficient estimator for ill-conditioned linear inverse problems
% with outliers', by Marta Martinez-Camara, Michael Muma, Abdelhak Zoubir
% and Martin Vetterli.
%
% Marta Martinez-Camara, LCAV-EPFL, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dependencies: the package CVX (http://cvxr.com/cvx/) is required.

%-----set up environment-----------------------------------
clear all; clc;
rng('shuffle'); % set the random seed to the clock (independent results each time this script runs)

% --- set up cvx ----------------------------------------------------
addpath('your path here'); % add your cvx path
cvx_setup; % set up cvx
cvx_quiet(true); % suppress output from cvx

%----set up parameters--------------------------------------------
m = 300;% Number of measurements
n = 120;% Number of unknowns
snr = 30;% SNR for the Gaussian noise
nReal = 100; % number of Monte Carlo realizations
out = [0 0.1 0.2 0.3 0.4]; % percentages of outliers to try
percOut = length(out);

nLmbds1 = 10; % how many lambdas for the LS and M estimator do we want to try
nLmbds2 = 7; % how many lambdas for the tau estimator do we want to try
lambda = zeros(5,nLmbds1);% Lambdas for LS loss function
lambdaH = zeros(5,nLmbds1);% Lambdas for M loss function
lambdaT = zeros(5,nLmbds2);% Lambdas for tau loss function

lambda(1,:) = logspace(-2,2,nLmbds1); % Range for the lambdas - out(1)
lambdaH(1,:) = logspace(-2,3,nLmbds1); 
lambdaT(1,:) = logspace(1,5,nLmbds2); 

lambda(2,:) = logspace(4,4.5,nLmbds1); % Range for the lambdas - out(2)
lambdaH(2,:) = logspace(0,2,nLmbds1);
lambdaT(2,:) = logspace(2,8,nLmbds2); 

lambda(3,:) = logspace(4.1,4.5,nLmbds1); % Range for the lambdas - out(3)
lambdaH(3,:) = logspace(1,3,nLmbds1); 
lambdaT (3,:)= logspace(3,8,nLmbds2); 

lambda(4,:) = logspace(4.3,4.5,nLmbds1); % Range for the lambdas - out(4)
lambdaH(4,:) = logspace(1,2,nLmbds1); 
lambdaT (4,:)= logspace(5,7,nLmbds2); 

lambda(5,:) = logspace(4.3,4.5,nLmbds1); % Range for the lambdas - out(5)
lambdaH(5,:) = logspace(1,2,nLmbds1); 
lambdaT (5,:)= logspace(6,8,nLmbds2); 

%%% Reg. tau estimator parameters
control.N  = 500;   % number of initial random solutions
control.k = 2;       % number of initial IRLS steps on each candidate 
control.t  = 5;      % number of best candidates to fully improve 
control.r = 2;       % number of iterations in scale approximation in case approx=1 
control.approx = 1 ; % if 0, fully compute S-scale, otherwise approximate

%---create synthetic data-------------------------------------
CondNumb = 1000; % condition number of the matrix
A = rand (m,n) ; % mixing matrix
x = zeros (n,1); 
x (30:42) = 1; % ground truth x
[U,S,V] = svd(A); % svd decomposition
S(S~=0)=linspace(CondNumb,1,min(m,n));
A = U*S*V'; % ill-conditioned matrix
y = A*x; % noisless measurements

%---storing the results------------------------------------
mseLS = zeros(percOut,nLmbds1,nReal);
mseT = zeros (percOut,nLmbds2,nReal);
mseH = zeros (percOut,nLmbds1,nReal);
mseTIK = zeros(nLmbds1,1);

for i = 1:percOut % percentage of outliers in the data
   display (['--------Outliers ', int2str(i), ' ----------'])
  for z = 1:nReal % number of realizations
    display (['Realization ', int2str(z), '.'])
    
    %--- Adding noise to the measurements-------------
    nout = round (out(i)*m);
    [dummysort,index]=sort(rand(m,1));
    index=index(1:nout);
    ynG = awgn (y,snr,'measured');% measurements with Gaussian noise
    yn = ynG;
    yn(index) = var(y)*10*rand (nout,1); % adding outliers


    %------Solving inverse problem---------------
    
    xLSs = zeros (n,nLmbds1);
    for k = 1:nLmbds1
       % LS loss function 
       xTIK = (A'*A + (lambda(i,k)^2)*eye(n))\A'*yn; % LS with Tik.
       mseTIK(k) = norm(xTIK - x);
       xLSs(:,k) = xTIK; 
       mseLS(i,k,z) = norm(xTIK - x);
       % M (Huber) loss function
       cvx_begin
        variable xxh(n);
        minimize(sum(huber((A*xxh-yn) /(1.483*mad(A*xTIK-yn,1)))) + lambdaH(i,k)*norm(xxh));% Huber loss function
       cvx_end
       xH = xxh;
       mseH(i,k,z) = norm(xH - x);
    end
    [opt, inx] = min(mseTIK);
    xLS = xLSs(:,inx);
    
    %tau loss function
    for k = 1:nLmbds2
     control.xHB = xLS; 
     result = FastTauRegFinal(A, yn,lambdaT(i,k),control); 
     xT = result.beta;
     mseT(i,k,z) = norm(xT - x);
    end
  end
end
%------Plotting results-------------------
minLSs = zeros(percOut,nReal);
minHs = zeros(percOut,nReal);
minTs = zeros(percOut,nReal);
for i = 1:percOut
minLSs(i,:) = min (squeeze(mseLS(i,:,:)),[],1);
minHs(i,:) = min (squeeze(mseH(i,:,:)),[],1);
minTs(i,:) = min (squeeze(mseT(i,:,:)),[],1);
end
avgMseLS = mean(minLSs,2);
avgMseH = mean(minHs,2);
avgMseT = mean(minTs,2);

plot(out,avgMseLS)
hold on
plot(out,avgMseH,'r')
plot(out,avgMseT,'g')

%---Saving resutls-----------------------
save('results.mat')
