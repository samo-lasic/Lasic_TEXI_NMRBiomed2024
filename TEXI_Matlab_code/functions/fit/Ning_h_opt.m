function opt = Ning_h_opt(opt)
% Makes sure that all needed fields in the options structure are present

opt.Ning_h.present = 1;

opt.Ning_h = msf_ensure_field(opt.Ning_h, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off'));
opt.Ning_h = msf_ensure_field(opt.Ning_h, 'do_plot', 0);
end