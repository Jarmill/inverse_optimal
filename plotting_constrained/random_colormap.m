function cmap = random_colormap(n, varargin)
%RANDOM_COLORMAP returns a random colormap.
%   cmap = RANDOM_COLORMAP(n)
%   cmap = RANDOM_COLORMAP(n, seed)
%
%   Inputs;
%   n ~ number of elements in the colormap
%   seed ~ seed color, around which cmap will be centered
%   
%   Outputs;
%   cmap ~ colormap (nx3)

% If seed is given
if length(varargin) >= 1
    
    % Get seed
    seed = varargin{1};
    if ~isequal(size(seed), [1, 3])
        error('seed must be a 1x3 color vector.');
    end
    if ~all(0 <= seed & seed <= 1)
        error('seed elements must be between 0 and 1');
    end
    
    % Minimum offsets along each component of seed
    min_off = min(seed, 1-seed);
    
    % Generate a random offset
    offset = min_off .* rand(1, 3);
    
    % For each component decide randomly whether to interpolate between the two
    % colors or not, by generating a binary string of 3 
    interpolate_random = randi([0, 1], 1, 3);

    % If all are zero, choose one randomly to set to 1
    if ~any(interpolate_random)
        interpolate_random( randi([1, 3], 1) ) = 1;
    end
    
    % For each component interpolate if corresponding interpolate_random is 1
    if interpolate_random(1)
        R = linspace(-offset(1), offset(1), n).' + seed(1);
    else
        R = seed(1) * ones(n, 1);
    end
    if interpolate_random(2)
        G = linspace(-offset(2), offset(2), n).' + seed(2);
    else
        G = seed(2) * ones(n, 1);
    end
    if interpolate_random(3)
        B = linspace(-offset(3), offset(3), n).' + seed(3);
    else
        B = seed(3) * ones(n, 1);
    end

    % Get the cmap
    cmap = [R, G, B];
    
    
% If seed is not given
else

    % Generate 2 random colors
    Ca = rand(1, 3);
    Cb = rand(1, 3);

    % Choose the concatenation of the two colors to be component-wise monotonic
    C1 = min(Ca, Cb);
    C2 = max(Ca, Cb);

    % For each component decide randomly whether to interpolate between the two
    % colors or not, by generating a binary string of 3 
    interpolate_random = randi([0, 1], 1, 3);

    % If all are zero, choose one randomly to set to 1
    if ~any(interpolate_random)
        interpolate_random( randi([1, 3], 1) ) = 1;
    end

    % For each component interpolate if corresponding interpolate_random is 1
    if interpolate_random(1)
        R = linspace(C1(1), C2(1), n).';
    else
%         R = (C1(1) + C2(1))/2 * ones(n, 1);
        R = rand(1, 1) * ones(n, 1);
    end
    if interpolate_random(2)
        G = linspace(C1(2), C2(2), n).';
    else
%         G = (C1(2) + C2(2))/2 * ones(n, 1);
        G = rand(1, 1) * ones(n, 1);
    end
    if interpolate_random(3)
        B = linspace(C1(3), C2(3), n).';
    else
%         B = (C1(3) + C2(3))/2 * ones(n, 1);
        B = rand(1, 1) * ones(n, 1);
    end

    % Get the cmap
    cmap = [R, G, B];
    
end

end