function blockMaster = gridGenerator(neurons, filterSize)

% determine grid dimensions
gridSize = sqrt(neurons);

% create matrix with grids
blockMaster = zeros(neurons, neurons);
c = 1;
x = zeros(gridSize, gridSize);
x(end - (filterSize - 1):end, end - (filterSize - 1):end) = 1;
x = circshift(x, [1, 1]);
for i = 1:gridSize 
    for j = 1:gridSize 
        temp = circshift(x, [i, j]);
        blockMaster(c,:) = temp(:)';
        c = c + 1; 
    end
end

% remove any redundant groups (for peace of mind!)
blockMaster = unique(blockMaster, 'rows', 'stable'); % remove redundant group vectors 
