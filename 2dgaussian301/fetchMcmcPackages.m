function fetchMcmcPackages()
    %Attempts to download dram package from 
    %http://www.helsinki.fi/~mjlaine/dram/dramcode.zip
    %and Spatial Econometrics toolbox from 
    %http://www.spatial-econometrics.com/html/jplv7.zip
    %Unzip it, and add to the path
    %These are required for MCMC-based error bars
    fprintf('Fetching dram package\n');
    thepath = fileparts(which('fetchMcmcPackages'));
    mkdir([thepath '/dram']);
    unzip('http://www.helsinki.fi/~mjlaine/dram/dramcode.zip',[thepath '/dram']);
    fprintf('Fetching spatial econometrics toolbox (6MB)\n');
    mkdir([thepath '/econo']);
    unzip('http://www.spatial-econometrics.com/html/jplv7.zip',[thepath '/econo']);
    
    addpath([thepath '/dram']);
    addpath([thepath '/dram/utils']);
    addpath(genpath([thepath '/econo']));
    fprintf('Fetched packages and added to path\n');
end