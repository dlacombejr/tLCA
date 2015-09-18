%Fit a Gaussian curve to data
xi = -10:10;
zi = exp(-(xi-3).^2/2/2^2) + .4*randn(size(xi));
clear opts

opts.positive = true; %So that the gain parameter is positive
opts.errorbars = 'none';
results = autoGaussianCurve(xi,zi,opts);

clf;
subplot(1,2,1);plot(xi(:),zi);
subplot(1,2,2);plot(xi(:),results.G);

%Look at the console to see the parameter estimates
results

%%
%All curve/surface fitting functions supports 2 types of error bars. 
%The first is based on bootstrapping; it's expensive but straightforward
%(ref: An Introduction to the Bootstrap, Efron & Tibshirani). 
%
%The second is based on MCMC (Markov-chain Monte Carlo). Specifically, the
%idea is to assume a model for the data (y = curve + normal noise) and perform 
%Bayesian inference on the model parameters. Because the model is 
%analytically intractable, inference is done by sampling from the
%posterior through MCMC. The MCMC method used is Delayed Rejection
%Adaptive Metropolis (DRAM) which is adaptive and deals quite well with the
%type of posterior involved here.
%
%MCMC is fast and accurate but it's more mathematically involved than
%bootstrapping. 


matlabpool open %call this once to enable parallel computing

%%

opts.errorbars = 'bootstrap';
tic;results = autoGaussianCurve(xi,zi,opts);toc;
results.quantiles %Prints quantiles of the points
                  %Look at results.quantiles.key for the corresponding
                  %probabilities. For example, 
                  %[results.quantiles.x0(1),results.quantiles.x0(end)] is a
                  %95% confidence interval for x0, while
                  %results.quantiles.sigmax(4) is the median of sigmax

%%
%MCMC based error bars require two packages to be installed
%Call this once to auto-download from the internet and add to path
fetchMcmcPackages();
%%
%Mcmc-based error bars
opts.errorbars = 'mcmc';
tic;results = autoGaussianCurve(xi,zi,opts);toc;
results.quantiles

clf;
subplot(2,2,1);plot(xi(:),zi);
subplot(2,2,2);plot(xi(:),results.G);
subplot(2,2,3);plot(xi(:),zi);
subplot(2,2,4);plot(xi(:),results.G);

%%
%Fit a Gaussian surface to data

[xi,yi] = meshgrid(-10:10,-20:20);
zi = exp(-(xi-3).^2/2/2^2-(yi+4).^2/2/3^2) + .2*randn(size(xi));
results = autoGaussianSurf(xi,yi,zi);

subplot(1,2,1);imagesc(xi(:),yi(:),zi);
subplot(1,2,2);imagesc(xi(:),yi(:),results.G);

%Look at the console to see the parameter estimates
results

%%
%Fit an isotropic Gaussian surface (sigmax == sigmay)
[xi,yi] = meshgrid(-10:10,-20:20);
zi = exp(-(xi-3).^2/2/2^2-(yi+4).^2/2/2^2) + .2*randn(size(xi));
clear opts;
opts.iso = true;
results = autoGaussianSurf(xi,yi,zi,opts);
results

subplot(1,2,1);imagesc(xi(:),yi(:),zi);
subplot(1,2,2);imagesc(xi(:),yi(:),results.G);

%%
%Fit a tilted Gaussian surface 
[xi,yi] = meshgrid(-10:10,-20:20);
xi = xi - 3;
yi = yi + 4;

angle = 45/180*pi;

xir = xi*cos(angle) + yi*sin(angle);
yir =-xi*sin(angle) + yi*cos(angle);
zi = exp(-(xir).^2/2/2^2-(yir).^2/2/3^2) + .2*randn(size(xi));

clear opts;
opts.tilted = true;
results = autoGaussianSurf(xi,yi,zi,opts);
results

subplot(1,2,1);imagesc(xi(:),yi(:),zi);
subplot(1,2,2);imagesc(xi(:),yi(:),results.G);


%%
%Fit a Gabor to the data
[xi,yi] = meshgrid(-0:1:20,-20:20);

opts.errorbars = 'mcmc';
xip = (xi-4)*cos(pi/6) + yi*sin(pi/6);
yip =-(xi-4)*sin(pi/6) + yi*cos(pi/6);
zi = exp(-(xip.^2+yip.^2)/2/2^2).*cos(xip*2*pi/5+pi/3) + .2*randn(size(xi));
results = autoGaborSurf(xi,yi,zi,opts);

subplot(1,2,1);imagesc(xi(:),yi(:),zi);
subplot(1,2,2);imagesc(xi(:),yi(:),results.G);

results
