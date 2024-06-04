function fn = make_seq_name(fn_prepend, seqName, fn_append)
c = {fn_prepend, seqName, fn_append};
c = c(~cellfun(@isempty,c));
fn = strjoin(c,'_');

