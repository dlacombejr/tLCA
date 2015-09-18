function [results] = autoGaborSurf(xi,yi,zi,opts)
    %function [results] = autoGaborSurf(xi,yi,zi,opts)
    %
    %Fit a Gabor (a sinusoid windowed by a Gaussian) to a surface. 
    %The Gabor is defined as follows:
    %zi = a*exp(-(xip.^2+yip.^2)/2/sigma^2)*cos(2*pi*xip/lambda + phase) + b
    %
    % Where:
    % xip = (xi-x0)*cos(theta) + (yi-i0)*sin(theta);
    % yip =-(xi-x0)*sin(theta) + (yi-i0)*cos(theta);
    %
    %On the assumption of Gaussian noise through maximum likelihood
    %(least-squares).
    %
    %The procedure is automatic in the sense that it is unlikely to get 
    %stuck in local minima which makes it appropriate in fully automated 
    %scenarios. This is accomplished by a trick: the absolute value of a
    %Gabor in the Fourier domain is a Gaussian. Thus the procedure finds
    %a good starting set of values for the frequency and orientation of the
    %Gabor by calling autoGaussianSurf on the Fourier transforms of the
    %data, uses that to get a good start guess for x0 and y0, and refines
    %the params using lsqcurvefit.
    %
    %Currently only regular grids (generated through meshgrid) are accepted
    %for xi and yi. 
    %
    %opts is a struct containing options:
    % opts.errorbars (none | bootstrap | mcmc) - Type of errors bars
    % requested (default:none, see TryAutoGaussianSurf for use)
    %
    %Example use:
    %
    % [xi,yi] = meshgrid(-10:10,-20:20);
    % xip = (xi-4)*cos(pi/6) + yi*sin(pi/6);
    % yip =-(xi-4)*sin(pi/6) + yi*cos(pi/6);
    % zi = exp(-(xip.^2+yip.^2)/2/2^2).*cos(xip*2*pi/5+pi/3) + .2*randn(size(xi));
    % results = autoGaborSurf(xi,yi,zi);
    %
    %See also autoGaussianSurf
    if nargin < 4
        opts = struct();
    end
    
    p = inputParser;
    p.KeepUnmatched = true;
    p.addOptional('errorbars','none');
    p.addOptional('bsIters',500);
    p.addOptional('mcmcIters',5000);
    p.addOptional('precision',.001*ones(8,1));
    p.parse(opts);
    opts = p.Results;
    
    %Start with fitting a Gaussian in the Fourier domain
    zif = fftshift(abs(fft2(zi)));
    %smooth a little, this helps
    g = [.2,1,.2];
    g = g'*g;
    zif = conv2(zif-mean(zif(:)),g,'same');
    
    ss = @(x) x(1:end-1);
    
    %Because of wraparound, do it in two shots
    [xip,yip] = meshgrid(linspace(-.5,0,size(xi,2)/2),ss(linspace(-.5,.5,size(xi,1)+1)));
    r1 = autoGaussianSurf(xip,yip,zif(:,1:size(xi,2)/2));
    
    %
    [xip,yip] = meshgrid(ss(linspace(-.5,.5,size(xi,2)+1)),linspace(-.5,0,size(xi,1)/2));
    r2 = autoGaussianSurf(xip,yip,zif(1:size(xi,1)/2,:));
    
    %Look at the quality of each fit to decide which to use
    
    if r1.sse < r2.sse
        %Go with r1
        coords = [r1.x0,r1.y0,max(r1.sigmax,r1.sigmay)];
    else
        coords = [r2.x0,r2.y0,max(r2.sigmax,r2.sigmay)];
    end
    
    %0 degrees is vertical
    theta = atan2(coords(2),coords(1));
    sf = min(sqrt(coords(1).^2+coords(2).^2),.5);
    dx = xi(1,2)-xi(1,1);
    dy = yi(2,1)-yi(1,1);
    da = .5*(dx+dy);
    lambda = 1/sf*da;
    sigma  = da/coords(3)/pi/2/1.3;
    phase = 0;
    
    %Get Gabor eval'd with current grid
    g = evalGabor([0,0,theta,lambda,sigma,phase],[xi(:)-mean(xi(:)),yi(:)-mean(yi(:))]);
    
    s = conv2(zi,reshape(g,size(xi)),'same');
    [y0,x0] = find(abs(s)==max(abs(s(:))));
    x0 = xi(1,x0);
    y0 = yi(y0,1);
    
    minsigma = da*.2;
    maxsigma = (max(xi(:))-min(xi(:))+max(yi(:))-min(yi(:)))/2*.4;
    
    %Now use lsqcurvefit to finish
    lb = [-Inf,-Inf,min(xi(:)),min(yi(:)),-Inf,-Inf,minsigma /1.01,-Inf]';
    ub = [ Inf, Inf,max(xi(:)),max(yi(:)), Inf, Inf,maxsigma + .01, Inf]';
    
    p0 = [x0,y0,theta,lambda,sigma,phase];
    g = evalGabor(p0,[xi(:),yi(:)]);
    
    p = [g,ones(size(g))]\zi(:);
    
    params0 = [p(1),p(2),p0]';
    varnames = {'a','b','x0','y0','theta','lambda','sigma','phase'};
    iscircular = [0,0,0,0,1,0,0,1]'==1;
    
    results = doFinalOptimization(@evalGabor,[xi(:),yi(:)],zi(:),params0,lb,ub,true(length(lb),1),varnames,iscircular,opts);
    
    %Collect the results
    results.G = reshape(results.G,size(zi));
end


function g = evalGabor(ps,X)
    xi = X(:,1);
    yi = X(:,2);
    x0 = ps(1);
    y0 = ps(2);
    theta = ps(3);
    lambda = ps(4);
    sigma = ps(5);
    phase = ps(6);
    
    xip =  (xi-x0)*cos(theta) + (yi-y0)*sin(theta);
    yip = -(xi-x0)*sin(theta) + (yi-y0)*cos(theta);
    
    g = exp(-(xip.^2+yip.^2)/2/sigma^2).*cos(2*pi*xip/lambda+phase);
end