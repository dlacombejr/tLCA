function [results] = autoGaussianCurve(xi,zi,opts)
    %function [results] = autoGaussianCurve(xi,zi,opts)
    %
    %Fit a curve zi = a*exp(-((xi-x0).^2/2/sigmax^2 + b
    %
    %On the assumption of Gaussian noise through maximum likelihood
    %(least-squares).
    %
    %The procedure is "automatic" in the sense that it is unlikely to get 
    %stuck in local minima which makes it appropriate in fully automated 
    %scenarios. This is accomplished through an initial exhaustive search 
    %for the parameters, followed by refinement with lsqcurvefit
    %
    %Currently only regular grids (generated through :) are accepted
    %for xi. 
    %Example use:
    %
    %xi = -10:10;
    %zi = exp(-(xi-3).^2) + randn(size(xi));
    %results = autoGaussianCurve(xi,zi)
    %
    %opts contains options:
    % opts.positive: if true, a must be > 0 (default false)
    % opts.bootstrap: if true, obtain error bars for parameters through
    %                 bootstrapping (default = false).
    %
    %Returns results, a struct with elements a,b,x0,sigmax for
    %the estimated model parameters, G which is the Gaussian evaluated with
    %the estimated model params, and sse, sse0 and r2, which are the sum of
    %squared error for the fitted Gaussian, the sum of squared error for a
    %model with only an offset, and the R^2 proportion of variance
    %accounted for.    
    %
    %opts is a struct containing options:
    % opts.errorbars (none | bootstrap | mcmc) - Type of errors bars
    % requested (default:none)
    % opts.positive: if true, gain (a) is required to be > 0
    %

    if nargin < 3
        opts = struct();
    end
    
    %parse inputs
    p = inputParser;
    p.KeepUnmatched = true;
    p.addOptional('positive',false);
    p.addOptional('errorbars','none');
    p.addOptional('bsIters',500);
    p.addOptional('mcmcIters',5000);
    p.addOptional('precision',[.001,.001,.001,.001]');
    p.parse(opts);
    opts = p.Results;
    
    sz = size(zi);
    
    %Verify that the grid is regular
    if any(any(abs(diff(xi,2,2)) >=1e3*eps))
        error('xi is not a regular grid');
    end
    
    if any(size(zi)~=size(xi))
        error('xi and zi are not the same size');
    end
    
    if numel(zi) ~= length(zi)
        error('zi must be 1-dimensional');
    end
    
    xi = xi(:);
    boundx = [min(xi),max(xi)];
    
    %Find a minimum sigma based on number of elements, range of x and y
    rgx = diff(boundx);
    minsigma = rgx/sz(2)/5;
    maxsigma = rgx*.3;
   
    sigmas = exp(log(minsigma):.3:log(maxsigma));
    
    rgx = -ceil(sz(2)/2)+1:sz(2)/2;
    
    res = zeros(length(sigmas),5);
    
    %Run through all the different values for sigma
    for ii = 1:length(sigmas)
        thefiltx = exp(-rgx.^2/2/sigmas(ii));
        %convolve zi with gaussian filters and find the maximum response
        %(matched filtering)
        zi2 = reflectconv(zi,thefiltx);
        if opts.positive
            [~,pos] = max(zi2(:));
        else
            [~,pos] = max(abs(zi2(:)));
        end
        x0 = xi(pos);
        
        %Determine the residual error for the optimal x for this sigma
        G = exp(-((xi-x0).^2)/2/sigmas(ii)^2);
        X = [G,ones(length(G),1)];
        ps = X\zi(:);
        res(ii,:) = [sum((zi(:) - X*ps).^2),ps(:)',x0,sigmas(ii)];
    end
    
    %Find sigma with corresponding least error
    [~,optsigma] = min(res(:,1));
    params0 = res(optsigma,2:end)';
    
    varnames = {'a','b','x0','sigmax'};
    lb = [-Inf,-Inf,xi(1),minsigma /1.01]';
    ub = [ Inf, Inf,xi(end),maxsigma + .01]';
    if opts.positive
        lb(1) = 0;
    end
    results = doFinalOptimization(@pointgaussian,xi(:),zi(:),params0,lb,ub,true(4,1),varnames,false(4,1),opts);
    results.G = reshape(results.G,size(xi));
end

%Convolution with reflected edge handling
function A = reflectconv(A,f)
    A = conv(A,f,'same');
end

function [a] = pointgaussian(x,xdata)
    a = exp(-.5*(xdata-x(1)).^2*(1./x(2).^2));
end