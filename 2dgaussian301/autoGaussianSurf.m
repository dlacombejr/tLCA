function [results] = autoGaussianSurf(xi,yi,zi,opts)
    %function [results] = autoGaussianSurf(xi,yi,zi,opts)
    %
    %Fit a surface zi = a*exp(-((xi-x0).^2/2/sigmax^2 + ...
    %                           (yi-y0).^2/2/sigmay^2)) + b
    %
    %On the assumption of Gaussian noise through maximum likelihood
    %(least-squares).
    %
    %The procedure is "robust" in the sense that it is unlikely to get 
    %stuck in local minima which makes it appropriate in fully automated 
    %scenarios. This is accomplished through an initial exhaustive search 
    %for the parameters, followed by refinement with lsqcurvefit
    %
    %Currently only regular grids (generated through meshgrid) are accepted
    %for xi and yi. 
    %Example use:
    %
    %[xi,yi] = meshgrid(-10:10,-20:20);
    %zi = exp(-(xi-3).^2-(yi+4).^2) + randn(size(xi));
    %results = autoGaussianSurf(xi,yi,zi)
    %
    %Returns results, a struct with elements a,b,x0,y0,sigmax,sigmay for
    %the estimated model parameters, G which is the Gaussian evaluated with
    %the estimated model params, and sse, sse0 and r2, which are the sum of
    %squared error for the fitted Gaussian, the sum of squared error for a
    %model with only an offset, and the R^2 proportion of variance
    %accounted for.
    %
    %opts is a struct containing options:
    % opts.errorbars (none | bootstrap | mcmc) - Type of errors bars
    % requested (default:none)
    % opts.iso (true|false) - If true, Gaussian is isotropic (sigmax ==
    % sigmay, default false)
    % opts.tilted (true|false) - If true, Gaussian is tilted by an angle,
    % default false). In that case the definition of the Gaussian becomes:
    % 
    % xip = (xi-x0)*cos(theta) + (yi-i0)*sin(theta);
    % yip =-(xi-x0)*sin(theta) + (yi-i0)*cos(theta);
    % zi = a*exp(-(xip.^2/2/sigmax^2 + ...
    %              yip.^2/2/sigmay^2)) + b
    %
    % opts.positive: if true, gain (a) is required to be > 0

    if nargin < 4
        opts = struct();
    end
    
    %parse inputs
    p = inputParser;
    p.KeepUnmatched = true;
    p.addOptional('positive',false);
    p.addOptional('errorbars','none');
    p.addOptional('bsIters',500);
    p.addOptional('mcmcIters',5000);
    p.addOptional('precision',[.001,.001,.001,.001,.001,.001,.001]');
    p.addOptional('iso',false);
    p.addOptional('tilted',false);
    p.parse(opts);
    opts = p.Results;

    sz = size(zi);
    
    %Verify that the grid is regular
    if any(any(abs(diff(xi,2,2)) >=1e3*eps)) || any(any(diff(yi,2,1) >= 1e3*eps))
        error('xi or yi is not a regular grid');
    end
    
    if any(size(zi)~=size(xi)) || any(size(zi)~=size(yi))
        error('xi, yi and zi are not the same size');
    end
    
    if opts.tilted && opts.iso
        error('A Gaussian cannot be both isotropic and tilted');
    end
    
    xi = xi(:);
    yi = yi(:);
    boundx = [min(xi),max(xi)];
    boundy = [min(yi),max(yi)];
    
    %Find a minimum sigma based on number of elements, range of x and y
    rgx = diff(boundx);
    minsigmax = rgx/sz(2)/5;
    maxsigmax = rgx/2;

    rgy = diff(boundy);
    minsigmay = rgy/sz(1)/5;
    maxsigmay = rgy/2;
    
    minsigma = min(minsigmax,minsigmay);
    maxsigma = max(maxsigmax,maxsigmay);
    sigmas = exp(log(minsigma):.3:log(maxsigma));
    
    rgx = [0:sz(2)/2,-ceil(sz(2)/2)+1:-1]';
    rgy = [0:sz(1)/2,-ceil(sz(1)/2)+1:-1]';
    
    res = zeros(length(sigmas),7);
    
    %Run through all the different values for sigma
    for ii = 1:length(sigmas)
        thefiltx = exp(-rgx.^2/2/sigmas(ii));
        thefilty = exp(-rgy.^2/2/sigmas(ii));
        %convolve zi with gaussian filters and find the maximum response
        %(matched filtering)
        zi2 = reflectconv(reflectconv(zi,thefilty)',thefiltx)';
        if opts.positive
            [~,pos] = max(zi2(:));
        else
            [~,pos] = max(abs(zi2(:)));
        end
        x0 = xi(pos);
        y0 = yi(pos);
        %[y0,x0] = ind2sub(sz,pos);
        
        %Determine the residual error for the optimal x, y for this sigma
        G = exp(-((xi-x0).^2+(yi-y0).^2)/2/sigmas(ii)^2);
        X = [G,ones(length(G),1)];
        ps = X\zi(:);
        res(ii,:) = [sum((zi(:) - X*ps).^2),ps(:)',x0,y0,sigmas(ii),sigmas(ii)];
    end
    
    %Find sigma with corresponding least error
    [~,optsigma] = min(res(:,1));
    
    %Fit the parameters again through lsqcurvefit
    if opts.iso
        lb = [-Inf,-Inf,boundx(1),boundy(1),minsigmax /1.01]';
        ub = [ Inf, Inf,boundx(2),boundy(2),maxsigmax + .01]';
        params0 = res(optsigma,2:end-1)';
        varnames = {'a','b','x0','y0','sigma'};
        iscircular = false(5,1);
        thefun = @pointisogaussian;
        opts.precision = opts.precision(1:5);
    elseif ~opts.tilted
        lb = [-Inf,-Inf,boundx(1),boundy(1),minsigmax /1.01,minsigmay /1.01]';
        ub = [ Inf, Inf,boundx(2),boundy(2),maxsigmax + .01,maxsigmay + .01]';
        params0 = res(optsigma,2:end)';
        varnames = {'a','b','x0','y0','sigmax','sigmay'};
        iscircular = false(6,1);
        thefun = @pointgaussian;
        opts.precision = opts.precision(1:6);
    else
        %Fit a Gaussian to the power spectrum of the thing
        theta = getInitialAngle(reshape(xi,sz),reshape(yi,sz),reshape(zi,sz),opts);
        
        %Because of wraparound, do it in two shots
        lb = [-Inf,-Inf,boundx(1),boundy(1),minsigmax /1.01,minsigmay /1.01,-Inf]';
        ub = [ Inf, Inf,boundx(2),boundy(2),maxsigmax + .01,maxsigmay + .01,Inf]';
        params0 = [res(optsigma,2:end-2)';res(optsigma,end-1:end)'.*[1.3,0.7]';theta];
        varnames = {'a','b','x0','y0','sigmax','sigmay','angle'};
        iscircular = [false(6,1);true];
        thefun = @pointtiltedgaussian;
    end
    
    if opts.positive
        lb(1) = 0;
    end
    
    results = doFinalOptimization(thefun,[xi(:),yi(:)],zi(:),params0,lb,ub,true(length(lb),1),varnames,iscircular,opts);
    
    %Collect the results
    results.G = reshape(results.G,size(zi));
end

function theta = getInitialAngle(xi,yi,zi,opts)
    zif = fftshift(abs(fft2(zi)));
    %smooth a little, this helps
    g = [.2,1,.2];
    g = g'*g;
    zif = conv2(zif-mean(zif(:)),g,'same');

    ss = @(x) x(1:end-1);
    
    opts.tilted = false;
    opts.errorbars = 'none';
    [xip,yip] = meshgrid(linspace(-.5,0,size(xi,2)/2),ss(linspace(-.5,.5,size(xi,1)+1)));
    r1 = autoGaussianSurf(xip,yip,zif(:,1:size(xi,2)/2),opts);

    %
    [xip,yip] = meshgrid(ss(linspace(-.5,.5,size(xi,2)+1)),linspace(-.5,0,size(xi,1)/2));
    r2 = autoGaussianSurf(xip,yip,zif(1:size(xi,1)/2,:),opts);

    %Look at the quality of each fit to decide which to use

    if r1.sse < r2.sse
        %Go with r1
        coords = [r1.x0,r1.y0,max(r1.sigmax,r1.sigmay)];
    else
        coords = [r2.x0,r2.y0,max(r2.sigmax,r2.sigmay)];
    end

    %0 degrees is vertical
    theta = atan2(coords(2),coords(1));
end


function [thef] = pointisogaussian(x,xdata)
    thef = exp(-.5*sum(bsxfun(@minus,xdata,[x(1),x(2)]).^2,2)*(1./x(3).^2));
end

function [thef] = pointtiltedgaussian(x,xdata)
    xdat = bsxfun(@minus,xdata,[x(1),x(2)]);
    xdat = xdat*[cos(x(5)),-sin(x(5)); sin(x(5)),cos(x(5))];
    thef = exp(-.5*(xdat.^2)*(1./[x(3);x(4)].^2));
end

function [thef] = pointgaussian(x,xdata)
    thef = exp(-.5*(bsxfun(@minus,xdata,[x(1),x(2)]).^2)*(1./[x(3);x(4)].^2));
end

%Convolution with reflected edge handling
function A = reflectconv(A,f)
    A = bsxfun(@times,fft([A(end:-1:1,:);A;A(end:-1:1,:)]),fft([f(1:floor(end/2));zeros(length(f)*2,1);f(floor(end/2)+1:end)]));
    A = ifft(A);
    A = A(length(f)+1:2*length(f),:);
end