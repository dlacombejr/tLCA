function results = doFinalOptimization(curvefun,X,z,p0,lb,ub,mask,varnames,iscircular,opts)
    %Function internally called by auto* to perform the final optimization
    %or simulation based on initial parameters. This is what calls
    %lsqcurvefit, does bootstrapping or MCMC
    switch opts.errorbars
        case 'none'
            %Maximum-likelihood
            params = getMLestimate(curvefun,'Iter',p0,X,z,lb,ub,mask);
        case 'bootstrap'
            if matlabpool('size') == 0
                warning('doFinalOptimization:noparallel',...
                        'Parallel computing disabled. Try calling matlabpool open before using bootstrap estimates, much faster');
            end
            I = ceil(rand(size(X,1),opts.bsIters)*size(X,1));
            parfor ii = 1:opts.bsIters
                params(:,ii) = getMLestimate(curvefun,'Off',p0,X(I(:,ii),:),z(I(:,ii)),lb,ub,mask);
            end
        case 'mcmc'
            if ~exist('dramrun.m','file')
                thepath = fileparts(which('doFinalOptimization'));
                addpath(genpath([thepath '/dram']));
                addpath(genpath([thepath '/econo']));
                if ~exist('dramrun.m','file')
                    error('Mcmc packages not installed. Run fetchMcmcPackages and try again');
                end
            end
            
            model.ssfun = @(p,z) sum((p(1)*curvefun(p(3:end),X)+p(2)-z).^2);
            model.priorfun = @(p,x) sum(opts.precision.*p(:).^2);
            
            G = curvefun(p0(3:end),X);
            sigma2 = mean((G-z).^2);
            
            simparams.par0 = p0;
            simparams.sigma2 = sigma2;
            simparams.n0 = 1;
            simparams.n  = numel(z);
            simparams.bounds = [lb';ub'];

            mcmcopts.nsimu = opts.mcmcIters;
            mcmcopts.qcov = eye(length(p0))*sigma2*100;
            mcmcopts.adaptint = 10;
            
            nmin = 5000;
            params = [];
            s2samples = [];
            while size(params,1) < nmin
                %Run dram
                [diagnostics,paramsn,s2samplesn] = dramrun(model,z,simparams,mcmcopts);
                params = [params;paramsn];
                s2samples = [s2samples;s2samplesn];
                
                %Run coda and look at convergence diagnostics
                r = coda(params);
                %Get the number of iterations coda says I should plus a
                %little padding
                nmin = max(cell2mat({r.n}));
                
                mcmcopts.nsimu = ceil(nmin-size(params,1)+500);
                mcmcopts.qcov = diagnostics.qcov;
                simparams.par0 = params(end,:)';
            end
            
            maxburn = max(cell2mat({r.nburn}));
            
            params = params(maxburn+1+500:end,:);
            s2samples = s2samples(maxburn+1+500:end,:);
            
            params = [params';sqrt(s2samples)'];
            varnames{end+1} = 'sigmanoise';
            results.diagnostics = diagnostics;
            iscircular = [iscircular;false];
        otherwise
            error('Unknown method of obtaining error bars: %s',opts.errorbars);
    end
    
    if size(params,2) == 1
        for ii = 1:length(varnames)
            results.(varnames{ii}) = params(ii);
        end
    else
        results = getSummaryStatistics(params',varnames,iscircular);
        results.samples = params;
        params = median(params');
        
        Gp = 0;
        for ii = 1:size(results.samples,2)
            p = results.samples(:,ii);
            Gp = Gp + curvefun(p(3:end),X)*p(1) + p(2);
        end
        results.Gp = Gp/size(results.samples,2);
    end
    results.G = curvefun(params(3:end),X)*params(1) + params(2);
    results.sse = sum((z-results.G).^2);
    results.sse0 = sum((z-mean(z)).^2);
    results.r2 = 1-results.sse/results.sse0;

end

function p = getMLestimate(curvefun,display,p0,X,z,lb,ub,mask)
    opts = optimset('Display',display,'Jacobian','on');
    params = lsqcurvefit(@(p,X) dodiff(curvefun,p,p0,X,mask),p0(mask),X,z,lb(mask),ub(mask),opts);
    p = p0;
    p(mask) = params;
end

function [eta,deta] = dodiff(curvefun,p,p0,x,mask)
%Do differentiation of curve/surf wrt parameters
    ps = p0;
    ps(mask) = p;
    p = ps;
    a = p(1);
    b = p(2);
    p = p(3:end);
    eta = curvefun(p,x);
    deta = zeros(length(eta),length(p)+2);
    delta = 1e-5;
    for ii = 1:length(p)
        dp = zeros(size(p));
        dp(ii) = delta;
        deta(:,ii+2) = (curvefun(p+dp,x)-eta)/delta;
    end
    deta(:,2) = 1;
    deta(:,1) = eta;
    deta = deta(:,mask);
end