%% formulate the constrained optimum sdp
Q = {[1 0; 0 1], [1 0; 0 1], [1 0; 0 1], [1 0; 0 1]};
f = {[0;0], [0;1], [1; 1], [0.5; 0]};

% Q = {[0.5 0; 0 2], [3 -1; -1 2], [1 0; 0 1], [1 0; 0 2]};

% y = [2; 0.5];
y = [1; -0.5];


EQ = false;
INEQ = true;

if EQ
%     Aeq=[1 0];
%     beq=[0.5];

    Aeq = [1 1];
    beq = [1];

else
    Aeq = [];
    beq = [];
end


if INEQ
    A=[1 0; 1 1];
    b=[0.6; 1.5];
else
    A=[];
    b=[];
end

%% get the program
[Fd, objd, indexer] = sdp_full_cons_lin_yalmip(y, Q, f, A, b, Aeq, beq)

%% solve the program
opts = sdpsettings('solver', 'mosek', 'savesolveroutput', true, 'savesolverinput', true);
sol = optimize(Fd, objd, opts);

% model = export(Fd, objd, opts);

if sol.problem == 0
    M_rec = value(indexer.M);
    M_rec_clean = M_rec;
    M_rec_clean(abs(M_rec_clean) <= 1e-6) = 0;
    e_rec = sort(eig(M_rec), 'descend');
    dist_rec = sqrt(value(indexer.dist));
end

%% Plot Figures
nf = length(f);
x_star = cell(1, nf);
for i = 1:nf
    x_star{i} = Q{i} \ (f{i});
end

figure(5)
clf
hold all
n_mesh = 500;
h_combined_sol = plot_mathcal_G(Q, x_star, n_mesh, 2);
h_combined_sol.LineStyle = 'none';
h_combined_sol.FaceAlpha = 0.5;
h_sol = plot_x_star(x_star, 150, winter(nf));
h_lev = plot_level_sets(Q, x_star, 1000, 1, winter(nf));
for ii = 1 : length(h_lev)
    h_lev(ii).HandleVisibility = 'Off';
end
scatter(y(1), y(2), 600, 'xk', 'LineWidth',3)