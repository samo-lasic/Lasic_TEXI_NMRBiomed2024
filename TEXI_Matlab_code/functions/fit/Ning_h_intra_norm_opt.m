function opt = Ning_h_intra_norm_opt(opt)
% Makes sure that all needed fields in the options structure are present

opt.Ning_h_intra_norm.present = 1;

opt.Ning_h_intra_norm = msf_ensure_field(opt.Ning_h_intra_norm, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off'));
opt.Ning_h_intra_norm = msf_ensure_field(opt.Ning_h_intra_norm, 'do_plot', 0);
end