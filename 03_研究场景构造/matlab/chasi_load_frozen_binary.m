function scene=chasi_load_frozen_binary(root,maxGroups)
% Native MATLAB reader for MATLAB column-major float64 binary export.
if nargin<2 || isempty(maxGroups),maxGroups=inf;end
m=jsondecode(fileread(fullfile(root,'manifest.json')));
assert(strcmp(m.modelId,'srs_joint_snapshot_v2'));
for i=1:numel(m.binaryFiles)
 p=fullfile(root,m.binaryFiles(i).name);
 assert(isfile(p) && strcmpi(chasi_sha256_file(p),m.binaryFiles(i).sha256), ...
  'Frozen input hash mismatch: %s',p);
end
c=m.config;G=min(numel(m.resources),maxGroups);r=struct([]);
for g=1:G
 row=m.resources(g);K=row.users;L=row.trps;Nr=row.rx;Nt=row.tx;rank=row.rank;
 prefix=sprintf('g%02d_',g);shape=[K,L,Nr,Nt];
 item=struct();item.H_true=readArray(root,[prefix 'H_true'],shape);
 item.H_hat=readArray(root,[prefix 'H_hat'],shape);item.H=item.H_hat;
 item.C=readArray(root,[prefix 'C'],[K,L,Nr,Nr]);item.C_stat=item.C;
 item.servingMask=logical(readArray(root,[prefix 'serve'],[L,K]));
 item.ueIds=readArray(root,[prefix 'ueIds'],[K,1]);
 item.initial=readArray(root,[prefix 'initial'],[L,K,Nt,rank]);
 item.resourceWeightNu=row.resourceWeightNu;
 item.rb_count=row.rbCount;item.powerWattPerTrp=c.powerWattPerTrp;
 source=fullfile(root,[prefix 'sourceCov.real.bin']);
 if isfile(source),item.sourceCov=readArray(root,[prefix 'sourceCov'],[K,L,Nr,Nr,4]);
 else,item.sourceCov=[];end
 if isempty(r),r=item;else,r(g)=item;end
end
scene=struct('manifest',m,'config',c,'resources',r, ...
 'input_root',root,'input_manifest_sha256',chasi_sha256_file(fullfile(root,'manifest.json')));
end

function a=readArray(root,name,shape)
p=fullfile(root,[name '.real.bin']);fid=fopen(p,'rb','ieee-le');assert(fid>=0);
cleaner=onCleanup(@()fclose(fid));a=fread(fid,prod(shape),'double=>double');
assert(numel(a)==prod(shape),'Binary shape mismatch: %s',name);clear cleaner
ip=fullfile(root,[name '.imag.bin']);
if isfile(ip)
 fid=fopen(ip,'rb','ieee-le');assert(fid>=0);cleaner=onCleanup(@()fclose(fid));
 b=fread(fid,prod(shape),'double=>double');assert(numel(b)==numel(a));
 a=complex(a,b);
end
a=reshape(a,shape);
end
