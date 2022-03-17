function [SDP_new] = mix_to_sdp(SDP_old)
%MIX_TO_SDP convert the mixed (linear + psd) sdp into a
%program with only psd blocks.
%   This conversion assumes an extremely standard form of the SDP, that the
%   linear constraints, where are linear 'l' variables are defined before the
%   PSD variables, and there are no free variables 'f'.

n_l = SDP_old.blk{1, 2};

ncon = length(SDP_old.b);

A_conv = cell(n_l, 1);
C_conv = cell(n_l, 1);
blk_conv = cell(n_l, 2);
for i = 1:n_l
    A_conv{i} = SDP_out.A{1}(i, :);
    C_conv{i} = 0;
    blk_conv{i, 1} = 's';
    blk_conv{i, 2} = 1;
end
% A_new = 

A_new = vertcat(A_conv, {SDP_old.A{2}});
C_new = vertcat(C_conv, {SDP_old.C{2}});
blk_new = vertcat(blk_conv, {SDP_old.blk{end, :}});

%package up the new SDP
SDP_new = struct;
SDP_new.A = A_new;
SDP_new.C = C_new;
SDP_new.blk = blk_new;
SDP_new.b = SDP_old.b;
SDP_new.ops = SDP_old.ops;

end

