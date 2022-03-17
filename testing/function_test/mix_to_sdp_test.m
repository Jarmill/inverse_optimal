%% The current STRIDE algorithm can only handle SDP blocks.
%test out a utility routine to convert the mixed (linear + psd) sdp into a
%program with only psd blocks.

%Future work: convert SOCP to 2x2 SDP blocks (less efficient when solving,
%but I can't go into the guts of SDPNAL to fix it).

load("sdp_opt_5_dual.mat", 'SDP_out')

n_l = SDP_out.blk{1, 2};

ncon = length(SDP_out.b);

A_conv = cell(n_l, 1);
C_conv = cell(n_l, 1);
blk_conv = cell(n_l, 2);
for i = 1:n_l
%     A_conv{i} = sparse(1, i, -1, 1, ncon);
    A_conv{i} = SDP_out.A{1}(i, :);
    C_conv{i} = 0;
    blk_conv{i, 1} = 's';
    blk_conv{i, 2} = 1;
end
% A_new = 

A_new = vertcat(A_conv, {SDP_out.A{2}});
C_new = vertcat(C_conv, {SDP_out.C{2}});
blk_new = vertcat(blk_conv, {SDP_out.blk{end, :}});

%package up the new SDP
SDP_mix = struct;
SDP_mix.A = A_new;
SDP_mix.C = C_new;
SDP_mix.blk = blk_new;
SDP_mix.b = SDP_out.b;
SDP_mix.ops = SDP_out.ops;

%% now test out the converted SDP and compare the outputs


[obj_new,X_new,y_new,Z_new,info_new,runhist_new] = sdpt3(SDP_mix.blk,SDP_mix.A,SDP_mix.C,SDP_mix.b); %,OPTIONS,X0,y0,Z0)

save(