function [sig_array, d_array, k_array, Nd, Nk] = sim2array(signals, sim_name)
sub_sig_struct = signals.(sim_name);
sub_seq_names = fields(sub_sig_struct);

sub_seq_names = sub_seq_names(contains(sub_seq_names,'_k'));

% get d and k arrays
N = numel(sub_seq_names);
for n_sub_seq_name = 1:N;
    sub_seq_name = sub_seq_names{n_sub_seq_name};
    d_array(n_sub_seq_name) = str2num(char(extractBetween(sub_seq_name,'d','_')));
    k_array(n_sub_seq_name) = str2num(char(extractAfter(sub_seq_name,'k')));
end

k_array = unique(k_array);
d_array = unique(d_array);
Nk = length(k_array);
Nd = length(d_array);

for nk = 1:Nk
    for nd = 1:Nd
        sig_array(nk, nd, :) = sub_sig_struct.(sprintf('d%d_k%d', d_array(nd), k_array(nk)));
    end
end

end