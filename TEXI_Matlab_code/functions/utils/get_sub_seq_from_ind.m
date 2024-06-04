function sub_seq = get_sub_seq_from_ind(seq, sub_ind)
sub_seq.n = length(sub_ind);

exclude_f = fields(sub_seq);
f = fields(seq);
f = f(~ismember(f,exclude_f));

for n = 1:numel(f)
    %if numel(seq.(f{n})) < sub_seq.n || contains(f{n},{'_array'})
    if numel(seq.(f{n})) == 1 || contains(f{n},{'_array'})
        sub_seq.(f{n}) = seq.(f{n});
    else
        sub_seq.(f{n}) = seq.(f{n})(sub_ind);
    end
end

% caution: g2_array might not be updated 
sub_seq.ntm = numel(unique(sub_seq.tm_ind));
sub_seq.tm_array = unique(sub_seq.tm);
sub_seq.nb = numel(unique(sub_seq.b_ind));
sub_seq.nc = numel(unique(sub_seq.gc_ind));
sub_seq.gc_array = unique(sub_seq.gc);
sub_seq.qc_array = unique(sub_seq.qc);

% re-index
[~, ~, ind] = unique(sub_seq.gc_ind); sub_seq.gc_ind = ind';
[~, ~, ind] = unique(sub_seq.b_ind); sub_seq.b_ind = ind';
[~, ~, ind] = unique(sub_seq.tm_ind); sub_seq.tm_ind = ind';
if isfield(sub_seq,'s_ind')
    [~, ~, ind] = unique(sub_seq.s_ind); sub_seq.s_ind = ind';
end

end