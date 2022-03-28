clear all;
close all;
clc;

% Define error functions
rmse = @(a,b) sqrt(sum((a-b).^2, 'all') ./ numel(a));
nrmse = @(a, b) sqrt(sum((a-b).^2 ./ a.^2, 'all') ./ numel(a));

% Load data
num_trial = 1;
data_filename = sprintf('ForceData_Optim_Format_Single_Trajectory_No%03d.mat', num_trial);
load(data_filename);


% Extract measured forces
f_hum = cell2mat(y.');

% Extract number of samples
N = size(f_hum, 2);
% Extract dimensionality of problem
n = size(f_hum, 1);
% Extract number of cost functions
nf = length(Q{1});

% Time vector
Time = linspace(0, 100, N);

% Prealocate lsq parameters
f_lsq = zeros(n, N);
alpha_lsq = zeros(nf, N);
d_lsq = zeros(1, N);
% Prealocate qcqp parameters
f_qcqp = zeros(n, N);
alpha_qcqp = zeros(nf, N);
d_qcqp = zeros(1, N);

% Optimoptions
qpopt = optimoptions(@quadprog, 'Display', 'None');

% For each sample
for ii = 1 : N
% for ii = 1 : 1
    
    % Extract inequality constraints that are bound type
    xl = -b{ii}(1:size(f_hum, 1));
    xu = b{ii}(size(f_hum, 1)+1:end);
    
    % Generate a random cost function parametrization
    alpha0 = rand(nf, 1);
    alpha0 = alpha0 ./ sum(alpha0);
    % Generate a feasible initial guess by solving the direct problem
    Qmix = 0;
    phimix = 0;
    for jj = 1 : nf
        Qmix = Qmix + alpha0(jj) * Q{ii}{jj};
        phimix = phimix+ alpha0(jj) * phi{ii}{jj};
    end
    [f0,~,~,~,mult0] = quadprog(Qmix, phimix, [], [], Aeq{ii}, beq{ii}, xl, xu, [], qpopt);
    % Extract multipliers
    lambda0 = mult0.eqlin;
    nu_minus0 = mult0.lower;
    nu_plus0 = mult0.upper;
    % Form initial solutions
    z0_lsq = [alpha0;lambda0;nu_minus0;nu_plus0];
    z0_qcqp = [f0;alpha0;lambda0;nu_minus0;nu_plus0];
    
    % Do the Keshavaraz method
    [f_lsq(:, ii), alpha_lsq(:, ii), d_lsq(ii)] = LSQ_IO_constrained(Q{ii}, phi{ii}, Aeq{ii}, beq{ii}, [], [], xl, xu, y{ii}, z0_lsq);

    % Do the QCQP method
    [f_qcqp(:, ii), alpha_qcqp(:, ii), d_qcqp(ii)] = QCQP_IO_constrained(Q{ii}, phi{ii}, Aeq{ii}, beq{ii}, [], [], xl, xu, y{ii}, z0_qcqp) ;

end

%%
figure ;
sgtitle({sprintf("RMSE(Human, KESH) = %.4e [N]; NRMSE(Human, KESH) = %.4e [%]", rmse(f_hum, f_lsq), nrmse(f_hum, f_lsq)); ...
         sprintf("RMSE(Human, MIND) = %.4e [N]; NRMSE(Human, MIND) = %.4e [%]", rmse(f_hum, f_qcqp), nrmse(f_hum, f_qcqp))}, ...
         'FontSize', 8);
[h_ax_hum, h_hum] = plot_vector_quantities(Time, f_hum, []) ;
[h_ax_lsq, h_lsq] = plot_vector_quantities(Time, f_lsq, []) ;
[h_ax_qcqp, h_qcqp] = plot_vector_quantities(Time, f_qcqp, []) ;

%%
figure;
subplot(2,1,1)
plot(Time, d_lsq);
ylim([0 2000])
subplot(2,1,2)
plot(Time, d_qcqp);
ylim([0 2000])