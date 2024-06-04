
function m = expm_data2fit(signal, xps, opt)

if ~isfield(opt,'n_inits')
    opt.n_inits = 1;
end

signal = double(signal);

% k  = m(1);
% f1 = m(2);
% D1 = m(3);
% D2 = m(4);

t_guess    = [5.000     0.5 10.00 100.0];
t_lb       = [0.0       0.0 5.00 50.0];
t_ub       = [20.000    1.0 20.0 700.0];

unit_to_SI = [1.0 1.0 1e-10 1e-10];

    function s = fun(t,varargin)
        m = t .* unit_to_SI;
        s = expm_fit2data(m, xps);
        %s = s(ind);
    end


if opt.n_inits > 0
    t = msf_fit(@fun, signal, t_lb, t_ub, opt.n_inits, opt.expm.lsq_opts);
else
    t = lsqcurvefit(@fun, t_guess, [], signal, ...
        t_lb, t_ub, opt.expm.lsq_opts);
end


m = t .* unit_to_SI;


if (opt.expm.do_plot)
    signal_fit = expm_fit2data(m, xps);
    figure(1),clf
    x = (1:numel(signal))';
    plot(x,signal,'.',x,signal_fit,'o');
    pause(0.05);
end

end
