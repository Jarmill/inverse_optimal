%use isosurface/contour commands to plot the region

load('opt_5_2d.mat')

n = 2;
m = 5;

% y_feas = [-1; 1];
% y_infeas = [2; 1];



func = @(y) alpha_lp(y, Q, x_star);

% s_feas = alpha_lp(y_feas);
% s_infeas = alpha_lp(y_infeas);

% s = func([y_feas, y_infeas]);

% Y = [-1 , 1; 
%       2,  1;
%       -2.4, 0.4]';


%% try plotting a surface/region using matlab's automatic routines

xrange=[-4, 3];
yrange=[-2.5,2.5];

Nx = 80;
Ny = 50;
[XX, YY] = meshgrid(linspace(xrange(1), xrange(2),Nx), ...
    linspace(yrange(1), yrange(2),Ny));
SS = zeros(Ny, Nx);
for i = 1:Nx
    for j = 1:Ny
        SS(j, i) = func([XX(j, i);YY(j, i)]);
    end
end

% SS = SS';

%% plot
figure(66)
clf
contour(XX, YY, SS, 1)

figure(67)
fcontour(@(x,y) func([x;y]), [xrange, yrange])

% s = func(Y)

function success = alpha_lp(y_in, Q, x_star)
    npts = size(y_in, 2);
    
    success = zeros(1, npts);
    for i = 1:npts
        [alpha_out, exitflag] = feas_lp_uncons(y_in(:, i),  Q, x_star);
        success(i) = (exitflag == 1);
    end
    
end