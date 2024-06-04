function m = Ning_h_intra_norm_data2fit(signal, xps, opt)
% m = [MD V k]

% MD  = m(1);
% V   = m(2);
% k   = m(3);
% S0  = 1
% Vintra = m(4)

if ~isfield(opt,'n_inits')
    opt.n_inits = 1;
end

signal = double(signal);

t_guess    = [1.000 1.000 0.00 1.000];
t_lb       = [0.001 0.001 0.00 0.001];
t_ub       = [4.000 4.000 40.0 4.000];


unit_to_SI = [1e-9 1e-18 1.0 1e-18];


    function s = fun(t,varargin)
        m = t .* unit_to_SI;
        s = Ning_h_intra_norm_fit2data(m, xps);
    end


if opt.n_inits > 0
    t = msf_fit(@fun, signal, t_lb, t_ub, opt.n_inits, opt.Ning_h_intra_norm.lsq_opts);
else
    t = lsqcurvefit(@fun, t_guess, [], signal, ...
    t_lb, t_ub, opt.Ning_h_intra_norm.lsq_opts);
end

m = t .* unit_to_SI;


if (opt.Ning_h_intra_norm.do_plot)
    signal_fit = Ning_h_intra_norm_fit2data(m, xps);
    x = (1:numel(signal))';
    figure(1),clf
    plot(x,signal,'.',x,signal_fit,'o');
    pause(0.05);
end

end