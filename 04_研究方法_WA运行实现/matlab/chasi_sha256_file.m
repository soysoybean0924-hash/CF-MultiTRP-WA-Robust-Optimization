function digest=chasi_sha256_file(path)
%CHASI_SHA256_FILE SHA-256 of file bytes on Windows and Linux, without Java.
assert(isfile(path));
path=char(path);
assert(isempty(regexp(path,'[\r\n]','once')), ...
    'File path cannot contain a newline');
if ispc
    % cmd.exe expands percent signs even inside quotes; reject shell syntax.
    assert(isempty(regexp(path,'["%&|<>^!]','once')), ...
        'File path contains Windows shell metacharacters');
    [status,output]=system(['certutil -hashfile "' path '" SHA256']);
    assert(status==0,'certutil failed to hash input (status %d): %s', ...
        status,output);
    matches=regexp(output,'(?m)^\s*([0-9a-fA-F]{64})\s*$','tokens');
else
    % Single quotes protect spaces, Chinese characters and shell syntax.
    quote=char(39);
    escaped=strrep(path,quote,[quote char(34) quote char(34) quote]);
    [status,output]=system(['sha256sum -- ' quote escaped quote]);
    assert(status==0,'sha256sum failed to hash input (status %d): %s', ...
        status,output);
    % GNU sha256sum prefixes a backslash when escaping a backslash in a name.
    if startsWith(output,'\'),output=output(2:end);end
    matches=regexp(output,'^([0-9a-fA-F]{64})[ *]','tokens');
end
assert(numel(matches)==1,'Unexpected SHA-256 command output');
digest=lower(matches{1}{1});
assert(~isempty(regexp(digest,'^[0-9a-f]{64}$','once')));
end
