function [Vw, b_array] = find_suitable_Vw_and_b(gc_array, qc_array, tm_array, b_array_guess, delta_c, delta_pi, inlude_DDE, g_lim, delta_lim)

% for VT-FEXI sequence we require
%  delta2 = 3*(b - qc.^2.*tmp_array)./(b*Vw - 2*qc.^2/delta_c);
%  delta = (- t_pi + sqrt(t_pi^2 + 4*delta2) )/2;

% start with a b_array_guess 
% define minimum b from crusher q and tm, min_bmin
% make array from min_bmin:min(b_array_guess)
% try sequences with different min_bmin from the table and Vw from a wide interval 
% the low limit of Vw is chose by
%     if bmin > 0
%         Vw_min = max(2*qc_array.^2/delta_c/bmin);
%     else
%         Vw_min = 2/delta_c./max(tm_array + 2*delta_c/3);
%     end
% check sequences for feasibility and get the gradient and delta ranges
% delta, g1, g2 are real, g<g_lim, delta_lim<delta<delta_lim

% if there are no feasible sequences, reduce the Vw range
% when feasible sequences are found for different bmin, find the minimum delta_max and min_max_g, save Vw, bmin
% reaper until Vw doesn't change too much

% find min b and max Vw with min delta
% b_low_lim prevents adjusting b too low

clc;
tm_array = [min(tm_array) max(tm_array)]; % just use min and max values
tmp_array = tm_array + 2*delta_c/3;

% define minimum b from crusher q and tm 
%3*(b - qc.^2.*tmp_array)./(b*Vw - 2*qc.^2/delta_c);
min_bmin = max(qc_array.^2).*max(tmp_array);

bmax = max(b_array_guess);

max_bmin = min(b_array_guess);

display(['bmin (crusher) = ' num2str(min_bmin*1e-6)])

Vw_factor = 1e4; %start with a wide range
Vw_steps = 1e3; 
b_steps = 1e2; 

if max_bmin < min_bmin
    display(['bmin is too small for the crusher !!!'])
    Vw = [];
    b_array = [];
    return
end

bmin_array = linspace(min_bmin,max_bmin,b_steps); % to find min b

N_max_iter = 10;
min_max_delta = nan*ones(N_max_iter,1);
min_max_g = min_max_delta;
Vw = min_max_delta;
bmin = min_max_delta;

Vw_factor_out = Vw_factor * 10;

for n_iter = 1:N_max_iter

      [par, isFeasible, min_delta, max_delta, min_g, max_g] = scan_seq_vs_b_Vw(Vw_factor, Vw_steps, bmin_array, bmax, ...
        gc_array, qc_array, delta_c, delta_pi, tm_array, inlude_DDE, delta_lim, g_lim);

    ind = find(isFeasible);


    if isempty(ind)
        Vw_factor = Vw_factor*.8;

    else
        min_max_delta(n_iter) = min(max_delta(ind));
        ind = find(max_delta == min_max_delta(n_iter));
        min_max_g(n_iter) = min(max_g(ind));
        ind = find(max_g == min_max_g(n_iter));
        ind = ind(1);

        Vw(n_iter) = par.Vw(ind);
        bmin(n_iter) = par.bmin(ind);


        if abs(Vw_factor_out-par.Vw_factor(ind)) < 1
            break
        else
            Vw_factor_out = par.Vw_factor(ind);
            Vw_factor = 3*Vw_factor_out; % zoom out a bit

        end

    end

end

% plot(min_max_delta_iter,'o')
% plot(Vw_iter,'o')
% plot(min_max_g_iter,'o')
% plot(bmin_iter,'o')

ind = find(~isnan(bmin));
if isempty(ind)
    display('no feasible sequence')
    Vw = [];
    b_array = [];

else
    ind = find(min(min_max_delta(ind)) == min_max_delta);
    ind = ind(1);
    min_max_delta = min_max_delta(ind);
    max_g = min_max_g(ind);
    Vw = Vw(ind);
    bmin = bmin(ind);


    display(sprintf('Vw_factor = %.2f reached in %d iterrations', Vw_factor, n_iter))

    b_array = b_array_guess;
    [~, min_ind] = min(b_array_guess);
    b_array(min_ind) = bmin;

        if ~isequal(b_array, b_array_guess)
        display('----- b_array has been modified! -----')
    end
end


function [par, isFeasible, min_delta, max_delta, min_g, max_g] = scan_seq_vs_b_Vw(Vw_factor, Vw_steps, bmin_array, bmax, gc_array, qc_array, delta_c, delta_pi, tm_array, inlude_DDE, delta_lim, g_lim)
Nb = length(bmin_array);
NVw = Vw_steps;

c = 1;
for nb = 1:Nb
    bmin = bmin_array(nb);
    b_array_tmp = [bmin bmax];

    if bmin > 0
        Vw_min = max(2*qc_array.^2/delta_c/bmin);
    else
        Vw_min = 2/delta_c./max(tm_array + 2*delta_c/3);
    end


    Vw_factor_array = linspace(1,Vw_factor,Vw_steps);
    Vw = Vw_min*Vw_factor_array;

    for nVw = 1:Vw_steps
        par.nb(c) = nb;
        par.nVw(c) = nVw;
        par.bmin(c) = bmin;
        par.bmax(c) = bmax;
        par.Vw(c) = Vw(nVw);
        par.Vw_factor(c) = Vw_factor_array(nVw);
        c = c+1;
    end
end

min_delta = nan*ones(c-1,1);
max_delta = min_delta;
min_g = min_delta;
max_g = min_delta;
isFeasible = min_delta;


parfor n = 1:c-1
%for n = 1:c-1
    seq = make_seq_Vw_gc(gc_array, delta_c, delta_pi, [par.bmin(n) par.bmax(n)], tm_array, par.Vw(n), inlude_DDE);

    min_delta(n) = min(seq.delta);
    max_delta(n) = max(seq.delta);

    min_g(n) = min([seq.g1 seq.g2]);
    max_g(n) = max([seq.g1 seq.g2]);

    isFeasible(n) = isempty(find(imag(seq.delta))) & isempty(find(imag(seq.g1))) & isempty(find(imag(seq.g2)));

    isFeasible(n) = isFeasible(n) & ...
        max_g(n) < g_lim  & ...
        min_delta(n) >= delta_lim(1) & ...
        max_delta(n) <= delta_lim(2);

end


