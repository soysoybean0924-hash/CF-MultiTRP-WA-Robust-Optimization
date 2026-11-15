function out=chasi_v6_channel_statistics(path,scene,group,tti,ueIds,previousTti)
% Conditional mean, diffuse receive covariance and adjacent-time cross C.
% No H_true is observed or used. Units match the noise-normalized trace.
if nargin<5,ueIds=1:scene.config.users;end
if nargin<6,previousTti=[];end
assert(isfile(path) && group>=1 && group<=scene.config.srsHops && ...
    tti>=0 && tti==fix(tti));
if ~isempty(previousTti)
 assert(previousTti>=0 && previousTti<tti && previousTti==fix(previousTti));
end
m=matfile(path);meta=m.meta;
assert(strcmp(meta.traceKind,'generated_parameters'));
U=numel(ueIds);L=scene.config.trps;Nr=scene.config.ueRxPorts;
Nt=scene.config.txPorts;
Q=complex(zeros(Nr,Nr,U,L));
meanH=complex(zeros(U,L,Nr,Nt,'single'));
cross=complex(zeros(Nr,Nr,U,L));
previousQ=complex(zeros(Nr,Nr,U,L));
previousMean=complex(zeros(U,L,Nr,Nt,'single'));
conditionTransitionCount=zeros(U,L);
previousConditionTransitionCount=nan(U,L);
crossEpochMatched=false(U,L);
persistent cachedPath cachedSize cachedDate cachedFields
file=dir(path);
if isempty(cachedPath) || ~strcmp(cachedPath,path)
 cachedFields=m.siteFields;
 cachedPath=path;cachedSize=file.bytes;cachedDate=file.datenum;
else
 assert(file.bytes==cachedSize && file.datenum==cachedDate, ...
     'Frozen channel statistics field changed during reconstruction');
end
fields=cachedFields;
for k=1:U
 u=ueIds(k);row=m.links(u,:);dropRow=m.dropDraw(u,:,:);
 for l=1:L
  site=ceil(l/3);drop=reshape(dropRow(1,site,:),1,2);
  [ray,rx,tx,id]=chasi_v6_tr38901_rebuild_link( ...
      row(l),fields{site},drop,meta,scene.trp_xy(l,:),tti);
  [currentQ,mean,rayNow]=chasi_v6_tr38901_ray_prior( ...
      ray,rx,tx,meta.carrierHz,meta.groupCenterOffsetHz(group));
  Q(:,:,k,l)=currentQ;
  meanH(k,l,:,:)=single(reshape(mean,1,1,Nr,Nt));
  conditionTransitionCount(k,l)=id.conditionTransitionCount;
  if ~isempty(previousTti)
   [priorRay,priorRx,priorTx,priorId]= ...
       chasi_v6_tr38901_rebuild_link(row(l),fields{site}, ...
       drop,meta,scene.trp_xy(l,:),previousTti);
   [oldQ,oldMean,rayPrevious]=chasi_v6_tr38901_ray_prior( ...
        priorRay,priorRx,priorTx,meta.carrierHz, ...
        meta.groupCenterOffsetHz(group));
   previousQ(:,:,k,l)=oldQ;
   previousMean(k,l,:,:)=single(reshape(oldMean,1,1,Nr,Nt));
   previousConditionTransitionCount(k,l)= ...
       priorId.conditionTransitionCount;
   crossEpochMatched(k,l)=priorId.isLos==id.isLos && ...
       priorId.conditionTransitionCount==id.conditionTransitionCount;
   if crossEpochMatched(k,l)
    assert(size(rayPrevious,3)==size(rayNow,3));
    small=fields{site}.condition{double(id.isLos)+1}.small;
    clusters=fields{site}.condition{double(id.isLos)+1}.clusters;
    displacement=id.positionMeters-priorId.positionMeters;
    for n=find(ray.clusterIndex>0).'
     for q=1:4
      phaseColumn=88*clusters+n+(q-1)*20*clusters;
      assert(phaseColumn<=size(small.kx,2) && ...
          ray.clusterIndex(n)==ceil(n/20));
      rho=mean(cos(small.kx(:,phaseColumn)*displacement(1)+ ...
          small.ky(:,phaseColumn)*displacement(2)));
      phaseDifference=ray.spatialPhasesRad(n,q)- ...
          priorRay.spatialPhasesRad(n,q);
      cross(:,:,k,l)=cross(:,:,k,l)+ ...
          chasi_v6_tr38901_phase_cross_kernel(rho)* ...
          exp(-1i*phaseDifference)* ...
          rayNow(:,:,n,q)*rayPrevious(:,:,n,q)'/Nt;
     end
    end
   end
  end
 end
end
out=struct('priorRx',Q,'meanH',meanH,'crossRx',cross, ...
    'previousPriorRx',previousQ,'previousMeanH',previousMean, ...
    'conditionTransitionCount',conditionTransitionCount, ...
    'previousConditionTransitionCount', ...
        previousConditionTransitionCount, ...
    'crossEpochMatched',crossEpochMatched, ...
    'group',group,'absoluteTti',tti,'previousTti',previousTti, ...
    'ueIds',double(ueIds(:)));
end
