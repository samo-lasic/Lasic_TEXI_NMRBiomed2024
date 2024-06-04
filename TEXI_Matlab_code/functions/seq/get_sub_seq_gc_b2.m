function [sub_ind, sub_seq] = get_sub_seq_gc_b2(seq, nc, b_min, b_max)
% b limits in s/mm^2
b_min = b_min * 1e6; % SI units
b_max = b_max * 1e6; % SI units

sub_ind = find(seq.gc_array(nc) == seq.gc & seq.b2 >= b_min & seq.b2 <= b_max);

sub_seq.n = length(sub_ind);
sub_seq.nc = 1;

exclude_f = fields(sub_seq);
f = fields(seq);
f = f(~ismember(f,exclude_f));

for n = 1:numel(f)
    if numel(seq.(f{n})) < sub_seq.n || contains(f{n},{'_array'})
        sub_seq.(f{n}) = seq.(f{n});
    else
        sub_seq.(f{n}) = seq.(f{n})(sub_ind);
    end
end

if isfield(sub_seq,'b_array')
    sub_seq.b_array = seq.b_array(seq.b_array < b_max);
    sub_seq.nb = length(sub_seq.b_array);
end
% sub_seq = seq;
% sub_seq.n = length(sub_ind);
% sub_seq.nc = 1;
% sub_seq.gc = seq.gc(sub_ind);
% sub_seq.gc_ind = seq.gc_ind(sub_ind);
% sub_seq.qc = seq.qc(sub_ind);
% sub_seq.b = seq.b(sub_ind);
% sub_seq.b_ind = seq.b_ind(sub_ind);
% sub_seq.b1 = seq.b1(sub_ind);
% sub_seq.b2 = seq.b2(sub_ind);
% sub_seq.tm = seq.tm(sub_ind);
% sub_seq.tm_ind = seq.tm_ind(sub_ind);
% sub_seq.delta = seq.delta(sub_ind);
% sub_seq.g = seq.g(sub_ind);
% sub_seq.gc_array = seq.gc_array(nc);
% %sub_seq. = seq.b_array(sub_ind);
% %sub_seq. = seq.tm_array;
% sub_seq.delta_array = squeeze(seq.delta_array(nc,:,:));
% sub_seq.g_array = squeeze(seq.g_array(nc,:,:));

end