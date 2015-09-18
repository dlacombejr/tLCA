function master = lateral_connection_generator(neurons)

% define grid size
dim = sqrt(neurons);

% create list of all pairwise x-y coordinates 
x = zeros(dim * dim, 2); 
c = 1; 
for i = 1:dim
    for j = 1:dim
        x(c, :) = [i, j]; 
        c = c + 1; 
    end
end

% create distance matrix of each cell from the center of the matrix
center_index = ceil(neurons / 2);
center = x(center_index, :); 
temp = zeros(dim, 1); 
for j = 1:size(x, 1)
    temp(j) = norm(center - x(j, :)); 
end
temp = reshape(temp, [dim, dim]); 

% shift the center of the matrix (zero distance) to the bottom right corner
temp = circshift(temp, [center_index - 1, center_index - 1]); 

% create master matrix 
master = zeros(neurons, neurons);
c = 1; 
for i = 1:dim
    for j = 1:dim
        new = circshift(temp, [j, i]); 
        master(c, :) = new(:)'; 
        c = c + 1; 
    end
end
