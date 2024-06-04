function opt = Ning_Gamma_opt(opt)
% Makes sure that all needed fields in the options structure are present

opt.Ning_Gamma.present = 1;

opt.Ning_Gamma = msf_ensure_field(opt.Ning_Gamma, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off'));
opt.Ning_Gamma = msf_ensure_field(opt.Ning_Gamma, 'do_plot', 0);
end