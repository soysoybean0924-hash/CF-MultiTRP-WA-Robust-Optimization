function digest=chasi_sha256_file(path)
%CHASI_SHA256_FILE Windows certutil SHA-256 for immutable input/config.
% The local MATLAB installation has no Java runtime. Restrict metacharacters
% before passing a quoted repository path to the system command.
assert(isfile(path));
path=char(path);
assert(isempty(regexp(path,'["%&|<>^\r\n]','once')));
[status,output]=system(['certutil -hashfile "' path '" SHA256']);
assert(status==0,'certutil failed to hash input');
matches=regexp(output,'(?m)^\s*([0-9a-fA-F]{64})\s*$','tokens');
assert(numel(matches)==1,'Unexpected certutil output');
digest=lower(matches{1}{1});
end
