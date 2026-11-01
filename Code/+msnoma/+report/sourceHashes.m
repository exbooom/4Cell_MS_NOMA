function hashes = sourceHashes(codeRoot)
%SOURCEHASHES Byte-level SHA256 of this implementation and its tests.
files=dir(fullfile(codeRoot,'**','*.m'));
hashes=struct('file',{},'sha256',{});
for i=1:numel(files)
    path=fullfile(files(i).folder,files(i).name);
    relative=extractAfter(path,[codeRoot filesep]);
    if startsWith(relative,['MS_NOMA' filesep]) || startsWith(relative,['__MACOSX' filesep])
        continue;
    end
    fid=fopen(path,'rb');
    assert(fid>0);
    cleanup=onCleanup(@() fclose(fid));
    bytes=fread(fid,Inf,'*uint8');
    clear cleanup
    engine=java.security.MessageDigest.getInstance('SHA-256');
    engine.update(bytes);
    digest=typecast(engine.digest(),'uint8');
    hashes(end+1)=struct('file',strrep(relative,filesep,'/'), ...
        'sha256',lower(reshape(dec2hex(digest,2)',1,[]))); %#ok<AGROW>
end
end
