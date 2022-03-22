clear all;
close all;
clc;

% Define polyhedron in 2D
n = 3;
% Define lower bounds
lb = zeros(3, 1);

% Define upper bounds
ub = ones(3, 1);

% No inequalities
Acons = [];
bcons = [];

% Define x1 + x2 + x3 <= 1
Acons = ones(1, 3);
bcons = [1];

% % Define x1 + x2 + x3 <= 1 & -x1 + x2 + x3 <= 0
% Acons = ones(2, 3);
% Acons(2, 1) = -1;
% bcons = [1; 0];

% No equalities
Aeqcons = [];
beqcons = [];
% % Define x1 = 0.4
% Aeqcons = [1, 0];
% beqcons = 0.4;

% % Define x1 + x2 == 1
% db = 0.0001;
% Acons1 = ones(1, 2);
% bcons1 = 1+db;
% Acons2 = -ones(1, 2);
% bcons2 = -(1-db);

figure;
hold on;
plot_feasible_set_constrained(Aeqcons, beqcons, Acons, bcons, lb, ub)
off = .25;
xlim([lb(1)-off, ub(1)+off])
ylim([lb(2)-off, ub(2)+off])
view(3)