
function opt = expm_opt(opt)
% function opt = expm_opt(opt)
% Makes sure that all needed fields in the options structure are present

opt.expm.present = 1;

opt.expm = msf_ensure_field(opt.expm, 'tmp', 1);
opt.expm = msf_ensure_field(opt.expm, 'lsq_opts', ...
    optimoptions('lsqcurvefit', 'display', 'off'));
opt.expm = msf_ensure_field(opt.expm, 'do_plot', 0);

end
