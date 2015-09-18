%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------%
%
% Machine Perception and Cognitive Robotics Laboratory
%
%     Center for Complex Systems and Brain Sciences
%
%              Florida Atlantic University
%
%------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------%
%
% Locally Competitive Algorithms Demonstration
% Using natural images data, see:
% Rozell, Christopher J., et al.
% "Sparse coding via thresholding and
% local competition in neural circuits."
% Neural computation 20.10 (2008): 2526-2563.
%
%------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear the workspace
clear all
close all
clc

% load in the images
load('IMAGES.mat')
I = IMAGES;

% set environment parameters
neurons = 121;      % number of neurons
patch_size = 256;   % patch size
batch_size = 1000;  % batch size
thresh = 0.1;       % LCA threshold 
h = .005;           % learning rate
blocksize = 3;      % neighborhood size
maxIter = 1000;     % maximum number of iterations

% create block matrix for inhibition weights
blockMaster = gridGenerator(neurons, blocksize); % call function to create this [groups X neurons]
groups = size(blockMaster, 1); % define number of groups (equals neurons with step == 1)

% create group inhibition weight matrix
G2 = blockMaster * blockMaster'; 
G2 = G2 ./ max(max(G2)); 
G2 = G2 - eye(neurons); 

% create lateral inhibition weight matrix
G1 = lateral_connection_generator(neurons);
G1 = 1 ./ G1; 
G1(G1 == inf) = 0; 

% randomly initialize wieghts
W = randn(patch_size, neurons);

% loop over training instances
for j = 1:maxIter 
    
    % normalize the weights
    W = W * diag(1 ./ sqrt(sum(W .^ 2, 1)));
    
    % create batch
    X = create_batch(I, patch_size, batch_size); % [patch_size X examples]
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get activation values for first (b1) and second level (b2)
    b1 = W' * X; % [neurons X examples]
    b2 = (blockMaster * sqrt(b1 .^ 2)) / blocksize; % [groups X examples]
    
    % LCA at layer 2
    u2 = zeros(groups,batch_size);
    for i = 1:5
        a2 = u2 .* (abs(u2) > thresh);
        u2 = 0.9 * u2 + 0.01 * (b2 - G2 * a2);        
    end
    a2=u2.*(abs(u2) > thresh); % [groups, batch_size]
    
    % project activations back down to first layer
    a1 = blockMaster' * a2; % [neurons X batch_size]
    a1 = a1 .* b1; % weight by first level activation

    % LCA at layer 1
    u1 = a1;
    for l =1:10
        a1=u1.*(abs(u1) > thresh);
        u1 = 0.9 * u1 + 0.01 * (b1 - G1*a1);
    end
    a1=u1.*(abs(u1) > thresh); % [groups, batch_size]
    
    % update the wieghts
    W = W + h * ((X - W * a1) * a1');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % visualize the learning procedure
    imagesc(filterplot(W))
    colormap(gray) 
    axis equal off
    drawnow()
    
    % display iteration
    disp(j)
    
end

% save the weights
save tLCA_weights.mat W
