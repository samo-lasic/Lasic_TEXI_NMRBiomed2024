function opt = Ning_g_opt(opt)
% Makes sure that all needed fields in the options structure are present

opt.Ning_g.present = 1;

opt.Ning_g = msf_ensure_field(opt.Ning_g, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off'));
opt.Ning_g = msf_ensure_field(opt.Ning_g, 'do_plot', 0);
end