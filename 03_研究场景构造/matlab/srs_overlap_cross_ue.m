function Ckj=srs_overlap_cross_ue(data,g,l,k,j)
%SRS_OVERLAP_CROSS_UE Reconstruct one conditional cross-UE covariance block.
% No true channel or realized estimation error is used. Does not allocate
% the full K-by-K-by-TRP-by-resource covariance tensor.
assert(ismember(data.modelId,{'srs_overlap_snapshot_v1','srs_overlap_jakes_v3'}));
K=size(data.pilotId,1); Nr=size(data.priorRx,1);
assert(all([k j]>=1 & [k j]<=K) && l>=1 && l<=size(data.priorRx,4));
d=data.resources(g).srs_diagnostic;
G=double(d.pilotGram); N=double(d.noiseCov); T=size(G,1);
rows=double(data.config.ulPortToDlRowsOneBased(:));
M=cell(K,1); S=N;
for u=1:K
 m=complex(zeros(T,Nr));
 for p=1:numel(rows)
  m(:,rows(p))=m(:,rows(p))+G(:,data.pilotId(u,p))*sqrt(data.pilotEnergy(u,p));
 end
 M{u}=m;
 R=double(data.priorRx(:,:,u,l));
 S=S+m*R*m';
end
Rk=double(data.priorRx(:,:,k,l)); Rj=double(data.priorRx(:,:,j,l));
a=double(d.a);
Ckj=-a^2*(Rk*M{k}'/S)*M{j}*Rj;
if k==j
 Ckj=Ckj+Rk;
end
end
