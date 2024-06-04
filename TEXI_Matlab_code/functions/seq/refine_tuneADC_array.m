function [tuneADC_array range_ADC_array, tuneADC] = refine_tuneADC_array(tuneADC_array, N, seq, Nint, ADC, tol, constraints, show_refine)
% includes gradient constraint

% n_tunable = find_tunable(seq, Nint, ADC, tuneADC_array, tol);
n_tunable = find_tunable(seq, Nint, ADC, tuneADC_array, tol, constraints);

if show_refine
    figure(2),clf
    plot(tuneADC_array,n_tunable,'.','MarkerSize',12)
end

[max_n_tunable, ind] = max(n_tunable);
tuneADC = tuneADC_array(ind(1));

ind1 = find(tuneADC_array < tuneADC_array(ind));
ind2 = find(tuneADC_array > tuneADC_array(ind));

len1 = length(ind1);
len2 = length(ind2);

if len1 > 1
    ind1 = ind1(round(len1 * .5));
end
if len2 > 1
    ind2 = ind2(round(len2 * .5));
end

tuneADC_array = linspace(tuneADC_array(ind1),tuneADC_array(ind2),N);

range_ADC_array = max(tuneADC_array)-min(tuneADC_array);


end