function [sub_ind, sub_seq] = get_sub_seq(seq, nc)
sub_ind = find(ismember(seq.gc,seq.gc_array(nc)));

sub_seq.n = length(sub_ind);
sub_seq.nc = 1;

exclude_f = fields(sub_seq);
f = fields(seq);
f = f(~ismember(f,exclude_f));

for n = 1:numel(f)
    %f{n}
    if numel(seq.(f{n})) < sub_seq.n || contains(f{n},{'_array'})
        sub_seq.(f{n}) = seq.(f{n});
    else
        sub_seq.(f{n}) = seq.(f{n})(sub_ind);
    end
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