function m = Ning_Gamma_data2fit(signal, xps, opt)
% m = [MD V k S0]

% MD  = m(1);
% V   = m(2);
% k   = m(3);
% S0  = m(4:end);

if ~isfield(opt,'n_inits')
    opt.n_inits = 1;
end

signal = double(signal);


for n = 1:max(xps.tm_ind)
    S0(n) = max(signal(xps.tm_ind == n));
end

t_guess    = [1.000 1.000 0.00 ones(size(S0))];
t_lb       = [0.001 0.001 0.00 zeros(size(S0))];
t_ub       = [4.000 4.000 40.0 2*ones(size(S0))];

unit_to_SI = [1e-9 1e-18 1.0 S0];


    function s = fun(t,varargin)
        m = t .* unit_to_SI;
        s = Ning_Gamma_fit2data(m, xps);
    end


if opt.n_inits > 0
    t = msf_fit(@fun, signal, t_lb, t_ub, opt.n_inits, opt.Ning_Gamma.lsq_opts);
else
    t = lsqcurvefit(@fun, t_guess, [], signal(ind), ...
    t_lb, t_ub, opt.Ning_Gamma.lsq_opts);
end


m = t .* unit_to_SI;


if (opt.Ning_Gamma.do_plot)
    signal_fit = Ning_Gamma_fit2data(m, xps);
    x = (1:numel(signal))';
    figure(1),clf
    plot(x,signal,'.',x,signal_fit,'o');
    pause(0.05);
end

end