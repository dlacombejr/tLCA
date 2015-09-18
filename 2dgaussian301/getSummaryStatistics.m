function results = getSummaryStatistics(samples,names,iscircular)
    p = [0.025,.05,.25,.5,.75,.95,.975];
    for ii = 1:size(samples,2)
        if iscircular(ii) > 0
            x = mean(cos(samples(:,ii)));
            y = mean(sin(samples(:,ii)));
            
            r = sqrt(y.^2+x.^2);
            results.means.(names{ii}) = atan2(mean(y),mean(x));
            results.stds.(names{ii})  = sqrt(-2*log(r));
        else
            results.means.(names{ii}) = mean(samples(:,ii));
            results.stds.(names{ii})  =  std(samples(:,ii));
            results.quantiles.(names{ii}) = quantile(samples(:,ii),p);
        end
    end
    
    results.quantiles.key = p;
end