load('test3_data.mat');

[dist_fmincon,alpha_con, exitflag, output] = fmincon_cons_search(y, Q, phi, Asdp, bsdp, Aeqsdp, beqsdp);