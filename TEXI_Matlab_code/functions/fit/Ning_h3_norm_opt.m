function opt = Ning_h3_norm_opt(opt)
% Makes sure that all needed fields in the options structure are present

opt.Ning_h3_norm.present = 1;

opt.Ning_h3_norm = msf_ensure_field(opt.Ning_h3_norm, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off'));
opt.Ning_h3_norm = msf_ensure_field(opt.Ning_h3_norm, 'do_plot', 0);
end