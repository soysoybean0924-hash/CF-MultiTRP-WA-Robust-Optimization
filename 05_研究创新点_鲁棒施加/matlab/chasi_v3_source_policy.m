function [resources,audit]=chasi_v3_source_policy(data,groups,policy)
%CHASI_V3_SOURCE_POLICY Experimental source/location policy for v3 input.
% Physical C is never changed here. policyCov is an additional design term.
% Names: target_residual, pilot_pollution, observation_noise, aging_innovation.
assert(ismember(data.modelId,{'srs_overlap_jakes_v3','srs_csi_joint_v4','srs_csi_joint_v5','srs_csi_joint_v6_formal'}));
names={'target_residual','pilot_pollution','srs_observation_noise', ...
    'aging_innovation','csi_rs_observation_noise'};
locations={'service','beam','power'};
assert(all(groups>=1 & groups<=numel(data.resources)));
sourceCount=size(policy.enabled,1);
assert(ismember(sourceCount,[4 5]) && size(policy.enabled,2)==3);
if sourceCount==4,names{3}='observation_noise_including_csi_rs';end
% Supported experimental placements: residual->beam, contamination->service,
% noise->beam, age->service/beam/power. Unsupported combinations fail closed.
allowed=logical([0 1 0;1 0 0;0 1 0;1 1 1;0 1 0]);
allowed=allowed(1:sourceCount,:);
assert(~any(policy.enabled(:) & ~allowed(:)),'Unsupported source/location pair');
assert(isequal(size(policy.strength),[sourceCount,3]));
assert(all(isfinite(policy.strength(:)) & policy.strength(:)>=0 & policy.strength(:)<=1));
assert(policy.powerReserveFraction>=0 && policy.powerReserveFraction<1);
if ~isfield(policy,'strengthMode'),policy.strengthMode='position';end
assert(ismember(policy.strengthMode,{'uniform','position','link_source'}));
if strcmp(policy.strengthMode,'link_source')
    assert(isfield(policy,'calibration') && ...
        isfield(policy.calibration,'gainErrorScale') && ...
        policy.calibration.gainErrorScale>0);
end
resources=data.resources(groups); audit=cell(1,numel(groups));
for ii=1:numel(resources)
    resources(ii).policyCov=[];
    resources(ii).powerWattPerTrp=double(data.config.powerWattPerTrp);
end
for gi=1:numel(groups)
    g=groups(gi); r=resources(gi); s=srs_overlap_source_statistics(data,g);
    assert(s.maxRelativeTotalMismatch<1e-5,'Source reconstruction failed');
    [K,L,Nr,~,~]=size(s.covariance);
    C=double(r.C); originalMask=r.servingMask;
    assert(isequal(size(C),[K,L,Nr,Nr]));
    if sourceCount==5
        assert(isfield(s,'fiveSourceCovariance') && ...
            size(s.fiveSourceCovariance,5)==5 && ...
            s.maxFiveSourceRelativeMismatch<1e-5, ...
            'Five-source policy requires separate CSI-RS covariance');
        sourceCov=s.fiveSourceCovariance;
        for k=1:K,for l=1:L
            physical=reshape(C(k,l,:,:),Nr,Nr);
            reconstructed=reshape(sum(sourceCov(k,l,:,:,:),5),Nr,Nr);
            if isfield(s,'crossSourceCovariance')
                reconstructed=reconstructed+reshape( ...
                    s.crossSourceCovariance(k,l,:,:),Nr,Nr);
            end
            assert(norm(physical-reconstructed,'fro')/ ...
                max(norm(physical,'fro'),eps)<1e-5, ...
                'Five-source covariance does not reconstruct physical C');
        end,end
    else
        sourceCov=s.covariance;
    end
    traces=zeros(K,L,sourceCount); gain=zeros(K,L);
    P=complex(zeros(size(C)));
    ageScale=min(1,max(0,double(r.age_ms)/ ...
        (double(data.config.revisitTti)*double(data.config.ttiMs))));
    for k=1:K, for l=1:L
        gain(k,l)=sum(abs(reshape(double(r.H(k,l,:,:)),[],1)).^2);
        if strcmp(data.modelId,'srs_csi_joint_v6_formal')
            gain(k,l)=gain(k,l)/size(r.H,4);
        end
        for q=1:sourceCount
            Q=reshape(sourceCov(k,l,:,:,q),Nr,Nr);
            traces(k,l,q)=max(0,real(trace(Q)));
        end
    end, end
    total=sum(traces,3); share=traces./max(total,eps);
    lambda=zeros(K,L,sourceCount,3);
    for q=1:sourceCount,for loc=1:3
        if ~policy.enabled(q,loc),continue;end
        if strcmp(policy.strengthMode,'uniform')
            base=policy.strength(1,1);
        else
            base=policy.strength(q,loc);
        end
        if strcmp(policy.strengthMode,'link_source')
            % Calibrated scale is fixed on development groups; all inputs
            % below are decision-time estimate and source statistics.
            risk=total./max(gain,policy.calibration.gainErrorScale);
            base=base*min(1,max(0,share(:,:,q).*risk));
        end
        if q==4,base=base*ageScale;end
        lambda(:,:,q,loc)=min(1,max(0,base));
    end,end
    for k=1:K,for l=1:L,for q=1:sourceCount
        if policy.enabled(q,2)
            Q=reshape(sourceCov(k,l,:,:,q),Nr,Nr);
            updated=reshape(P(k,l,:,:),Nr,Nr)+lambda(k,l,q,2)*Q;
            P(k,l,:,:)=reshape(updated,1,1,Nr,Nr);
        end
    end,end,end
    r.policyCov=P;
    % One replacement per UE, preserving the number of serving TRPs.
    if any(policy.enabled(:,1))
        score=gain./max(gain+total,eps);
        for q=1:sourceCount
            if policy.enabled(q,1)
                score=score-lambda(:,:,q,1).*share(:,:,q);
            end
        end
        for k=1:K
            present=find(originalMask(:,k)); absent=find(~originalMask(:,k));
            if isempty(present)||isempty(absent),continue;end
            [~,i]=min(score(k,present)); [~,j]=max(score(k,absent));
            if score(k,absent(j))>score(k,present(i))+policy.serviceMargin
                r.servingMask(present(i),k)=false;
                r.servingMask(absent(j),k)=true;
            end
        end
    end
    cap=double(data.config.powerWattPerTrp);
    reserve=zeros(K,L);
    for q=1:sourceCount
        if policy.enabled(q,3)
            reserve=reserve+lambda(:,:,q,3).*share(:,:,q);
        end
    end
    % The existing solver exposes a per-TRP cap, not per-link power bounds.
    % Use worst served-link risk as an explicit conservative TRP reserve.
    served=double(r.servingMask');
    perTrp=max(reserve.*served,[],1);
    r.powerWattPerTrp=cap*(1-policy.powerReserveFraction*max(perTrp));
    resources(gi)=r;
    audit{gi}=struct('group',g,'sources',{names(1:sourceCount)},'locations',{locations}, ...
        'physicalCUnchanged',isequaln(double(r.C),C), ...
        'sourceMismatch',s.maxRelativeTotalMismatch, ...
        'ageMs',double(r.age_ms),'ageScale',ageScale, ...
        'sourceTrace',traces,'sourceShare',share, ...
        'sourceCrossCovariance',s.crossSourceCovariance, ...
        'appliedStrength',lambda, ...
        'serviceMaskBefore',originalMask,'serviceMaskAfter',r.servingMask, ...
        'policyCovariance',P,'powerReserveByTrp',perTrp, ...
        'physicalPowerCapW',cap,'effectivePowerCapW',r.powerWattPerTrp);
end
audit=[audit{:}];
end
