function m = determine_simplex_grid_partition(n,Ngrid)
%DETERMINE_SIMPLEX_GRID_PARTITION determines the largest possible number of
%grid points along each edge of a simplex, given the maximum total number
%of grid points allowed and the affine dimension of the simplex.
%
%   m = determine_simplex_grid_partition(n,Ngrid)
%   Inputs:
%   n ~ affine dimension of simplex
%   Ngrid ~ maximum number of total points on the simplex

%   Outputs:
%   m ~ number of grid points along a single edge of a simplex

% Lower and upper bounds based on bounds of the nchoosek function
lb = floor(nthroot(Ngrid * factorial(n), n)) + 1 - n;
ub = ceil(nthroot(Ngrid * factorial(n), n));

m = lb;
for mm = lb+1 : 1 : ub
    
    % Calculate total number of gridpoints
    Ntot = nchoosek(n+mm-1, n);
    
    % If total number of grid points is bigger than requested
    if Ntot > Ngrid
        break
    end
    
    % Else increment m
    m = mm;
end
end