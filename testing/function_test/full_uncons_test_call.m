load('opt_3_2d.mat')

y_out = [2; 1];

[SDP_out,recoverdata,diagnostic,interfacedata] = sdp_full_uncons(y_out, Q, x_star);