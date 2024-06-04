function s = msf_ensure_field(s, f, v)
% function s = msf_ensure_field(s, f, v)
% function from https://github.com/markus-nilsson/md-dmri

% Sets the field 'f' in the structure 's' to the value 'v', unless there is
% already a field 'f' in 's'.

if (~isfield(s, f)), s.(f) = v; end
end