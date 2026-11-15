function stats=srs_overlap_source_statistics(data,groupIndex)
%SRS_OVERLAP_SOURCE_STATISTICS Exact linear-posterior error attribution.
% For y=sum_j M_j x_j+n and A_k=Q_k,: M' / S, the estimation error is
% (I-A_k M_k)x_k - sum_{j~=k} A_k M_j x_j - A_k n.
% Target and pollution need a separate signed cross term when Q_k,j~=0.
% Aging adds independent innovation and scales measurement terms by a^2.
% This is conditional on the project's Gaussian prior/finite-pilot model;
% it is not a fit of Arvidson Eq.23 or an externally validated split.
assert(ismember(data.modelId,{'srs_overlap_jakes_v3','srs_csi_joint_v4','srs_csi_joint_v5','srs_csi_joint_v6_formal'}));
if ismember(data.modelId,{'srs_csi_joint_v4','srs_csi_joint_v5','srs_csi_joint_v6_formal'})
    r=data.resources(groupIndex);parts=double(r.sourceCov);
    [K,L,Nr,~,~]=size(parts);mismatch=zeros(K,L);
    for k=1:K,for l=1:L
        total=reshape(sum(parts(k,l,:,:,:),5),Nr,Nr);
        C=reshape(double(r.C(k,l,:,:)),Nr,Nr);
        if isfield(r,'sourceCross')
            total=total+reshape(double(r.sourceCross(k,l,:,:)),Nr,Nr);
        end
        mismatch(k,l)=norm(total-C,'fro')/max(norm(C,'fro'),eps);
    end,end
    stats=struct('group',groupIndex,'ueIds',double(r.ueIds(:)), ...
        'componentNames',{{'target_residual','pilot_pollution', ...
        'observation_noise_including_csi_rs','aging_innovation'}}, ...
        'covariance',parts,'maxRelativeTotalMismatch',max(mismatch,[],'all'), ...
        'relativeTotalMismatch',mismatch,'scope', ...
        'exact stored linear Gaussian sequential-update attribution', ...
        'ageMs',double(r.age_ms));
    stats.crossSourceCovariance=complex(zeros(K,L,Nr,Nr));
    if isfield(r,'sourceCross')
        stats.crossSourceCovariance=double(r.sourceCross);
    end
    if isfield(r,'sourceCov5')
        five=double(r.sourceCov5);fiveMismatch=zeros(K,L);
        for k=1:K,for l=1:L
            total=reshape(sum(five(k,l,:,:,:),5),Nr,Nr);
            C=reshape(double(r.C(k,l,:,:)),Nr,Nr);
            if isfield(r,'sourceCross')
                total=total+reshape(double(r.sourceCross(k,l,:,:)),Nr,Nr);
            end
            fiveMismatch(k,l)=norm(total-C,'fro')/max(norm(C,'fro'),eps);
        end,end
        stats.fiveSourceNames={'target_residual','pilot_pollution', ...
            'srs_observation_noise','aging_innovation', ...
            'csi_rs_observation_noise'};
        stats.fiveSourceCovariance=five;
        stats.fiveSourceRelativeMismatch=fiveMismatch;
        stats.maxFiveSourceRelativeMismatch=max(fiveMismatch,[],'all');
        stats.conditionalErrorMean=complex(zeros(K,L,Nr, ...
            size(r.H,4)));
        stats.crossSourceCovarianceAssumption= ...
            'target and pollution share cross-UE prior; stored signed cross plus five PSD terms reconstruct physical C';
        maximumCross=0;
        for k=1:K,for l=1:L
            maximumCross=max(maximumCross,norm(reshape( ...
                stats.crossSourceCovariance(k,l,:,:),Nr,Nr),'fro'));
        end,end
        stats.crossSourceCovarianceMaxFro=maximumCross;
    end
    return
end
r=data.resources(groupIndex); pilot=double(data.pilotId);
energy=double(data.pilotEnergy); prior=double(data.priorRx);
jointPrior=[];
if isfield(data,'jointPrior') && ~isempty(data.jointPrior)
    jointPrior=double(data.jointPrior);
end
gram=double(r.srs_diagnostic.pilotGram);
noiseCov=double(data.pilotNoiseVariance)*gram;
rows=double(r.srs_rows); a=double(r.srs_diagnostic.a);
[K,L,Nr,~]=size(r.H); allK=size(pilot,1); T=size(gram,1);
assert(numel(rows)==size(pilot,2) && size(prior,3)==allK);
component=complex(zeros(K,L,Nr,Nr,4));
crossComponent=complex(zeros(K,L,Nr,Nr));
relativeMismatch=zeros(K,L);
M=complex(zeros(T,Nr,allK));
for j=1:allK
    for p=1:numel(rows)
        M(:,rows(p),j)=M(:,rows(p),j)+ ...
            sqrt(energy(j,p))*gram(:,pilot(j,p));
    end
end
for l=1:L
    if ~isempty(jointPrior)
        assert(size(jointPrior,1)==allK*Nr && ...
            size(jointPrior,2)==allK*Nr && ...
            size(jointPrior,3)==L);
        Q=jointPrior(:,:,l);Q=(Q+Q')/2;
    else
        Q=complex(zeros(allK*Nr));
        for j=1:allK
            ix=(j-1)*Nr+(1:Nr);
            Q(ix,ix)=prior(:,:,j,l);
        end
    end
    measure=complex(zeros(T));
    fullM=reshape(M,T,allK*Nr);
    measure=fullM*Q*fullM';
    S=(measure+noiseCov); S=(S+S')/2;
    for k=1:K
        id=double(r.ueIds(k)); J=M(:,:,id); R=prior(:,:,id,l);
        ix=(id-1)*Nr+(1:Nr);
        other=setdiff(1:allK*Nr,ix);
        A=(Q(ix,:)*fullM')/S;
        own=eye(Nr)-A*J;
        ownCov=own*R*own';
        pollution=A*fullM(:,other)*Q(other,other)* ...
            fullM(:,other)'*A';
        cross=-own*Q(ix,other)*fullM(:,other)'*A'- ...
            A*fullM(:,other)*Q(other,ix)*own';
        observation=A*noiseCov*A';
        innovation=(1-a^2)*R;
        parts=cat(3,a^2*ownCov,a^2*pollution, ...
            a^2*observation,innovation);
        for q=1:4
            C=(parts(:,:,q)+parts(:,:,q)')/2;
            component(k,l,:,:,q)=reshape(C,1,1,Nr,Nr,1);
        end
        crossComponent(k,l,:,:)=reshape(a^2*(cross+cross')/2, ...
            1,1,Nr,Nr);
        total=sum(parts,3)+a^2*cross;
        stored=reshape(double(r.C(k,l,:,:)),Nr,Nr);
        relativeMismatch(k,l)=norm(total-stored,'fro')/ ...
            max(norm(stored,'fro'),eps);
    end
end
stats=struct('group',groupIndex,'ueIds',double(r.ueIds(:)), ...
    'componentNames',{{'target_residual','pilot_pollution', ...
        'observation_noise','aging_innovation'}}, ...
    'covariance',component,'maxRelativeTotalMismatch', ...
        max(relativeMismatch,[],'all'), ...
    'crossSourceCovariance',crossComponent, ...
    'relativeTotalMismatch',relativeMismatch, ...
    'scope','exact for assumed linear Gaussian model; not external validation', ...
    'ageMs',double(r.age_ms));
end
