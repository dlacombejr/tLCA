% add gabor fit library to general path
path = '2dgaussian301/'; 
addpath(genpath(path))

% load in the wieghts
load tLCA_weights.mat W
W = W';

% determine number of neurons and dimensions (assuming square)
[neurons, dim_squared] = size(W); 
dim = sqrt(dim_squared); 

% set up meshgrid
[xi,yi] = meshgrid(1:dim, 1:dim);

% fit gabors and save in master matrix
fits = zeros(neurons, dim_squared); 
freq = zeros(neurons, 1);
orient = zeros(neurons, 1);
phase = zeros(neurons, 1);

for i = 1:neurons
    
    % grab current receptive field
    zi = reshape(W(i, :), [dim, dim]); 
    
    % normalization
    zi = zi ./ max(max(abs(zi))); 
    
    % get results
    results = autoGaborSurf(xi,yi,zi);
    
%     % plot findings
%     subplot(1, 2, 1)
%     imagesc(zi)
%     subplot(1, 2, 2)
%     imagesc(results.G); 
%     pause
        
    % append results into vectors
    fits(i, :) = results.G(:)'; 
    freq(i) = results.lambda; 
    orient(i) = results.theta; 
    phase(i) = results.phase; 
    
end

% put the results in a structure
f1 = 'fits';
f2 = 'freq';
f3 = 'orient';
f4 = 'phase';
results = struct(f1, fits, f2, freq, f3, orient, f4, phase); 
