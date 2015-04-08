%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to generate Figure 3 from the paper 'A new robust 
% and efficient estimator for ill-conditioned linear inverse problems
% with outliers', by Marta Martinez-Camara, Michael Muma, Abdelhak Zoubir
% and Martin Vetterli.
%
% You need cvx to run it (http://cvxr.com/cvx/). Set up your cvx path before using it.
%
% Marta Martinez-Camara, LCAV-EPFL, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------Set up environment------------
clear all; clc;
rng('shuffle');
addpath('your path here'); % add the cvx path
cvx_setup; % set up cvx
cvx_quiet(true); % suppress output from cvx
% Generating synthetic data to test the algorithm...
m = 300;
n = 120;
snr = 30;
nReal = 200;

%------Generating data------------
A = 10*rand (m,n) ; % synthetic matrix
x = zeros (n,1); 
x (30:42) = 1; % syntehtic source
y = A*x; % noisless measurements
ynG = awgn (y,snr,'measured');
out = [0 0.1 0.2 0.3 0.4 ]; % percentage of out to try
percOut = length(out);
mseLS = zeros(percOut,nReal);
mseT = zeros (percOut,nReal);
mseH = zeros (percOut,nReal);

for i = 1:percOut; % percentage of outliers in the data
  display (['--------Outliers ', int2str(i), ' ----------'])
  for z = 1:nReal % number of realizations
    display (['Realization ', int2str(z), '.'])
  nout = round (out(i)*m);
  [dummysort,index]=sort(rand(m,1));
  index=index(1:nout);
  yn = ynG;
  
  yn(index) = var(y)*10*rand (nout,1); % several ones


  %------Solving inverse problem---------------
  % least squares solution
  xLS = pinv(A)*yn; 
  control.N  = 500;
  control.k = 2;
  control.t  = 5;
  control.r = 2;
  control.approx = 1 ;
  % tau estimate
  result = FastTauFinal(A,yn,control);
  % M (Huber) estimate
  cvx_begin
    variable xh(n);
    minimize(sum(huber((A*xh-yn)/(1.483*mad(A*xLS-yn,1)))));% Huber loss function
  cvx_end
  xH = xh;
  % storing results
  mseLS(i,z) = norm(xLS - x);
  mseT(i,z) = norm(result.beta -x);
  mseH(i,z) = norm(xH -x);
  end
end

avgMseLS = mean(mseLS,2);
avgMseH = mean(mseH,2);
avgMseT = mean(mseT,2);

%------Plot solutions-------------------
plot(out,avgMseLS)
hold on
plot(out,avgMseH,'r')
plot(out,avgMseT,'g')

%------Save results---------------
save('results.mat')
