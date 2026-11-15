function cal=chasi_v3_calibrate_strength(data,developmentGroups)
% Lock strength normalization using development groups only.
assert(ismember(data.modelId,{'srs_overlap_jakes_v3','srs_csi_joint_v4','srs_csi_joint_v5','srs_csi_joint_v6_formal'}));
assert(~isempty(developmentGroups) && all(developmentGroups>=1) && ...
    all(developmentGroups<=numel(data.resources)));
gains=[];ratios=[];
for g=developmentGroups(:)'
    r=data.resources(g);[K,L,Nr,Nt]=size(r.H);
    for k=1:K,for l=1:L
        h=reshape(double(r.H(k,l,:,:)),[],1);
        gain=sum(abs(h).^2);
        if strcmp(data.modelId,'srs_csi_joint_v6_formal')
            gain=gain/Nt; % Match the per-TX-column RX covariance unit.
        end
        gains(end+1)=gain; %#ok<AGROW>
        if ~strcmp(data.modelId,'srs_csi_joint_v6_formal')
            e=reshape(double(r.H_true(k,l,:,:))-double(r.H(k,l,:,:)),[],1);
            C=reshape(double(r.C(k,l,:,:)),Nr,Nr);
            ratios(end+1)=mean(abs(e).^2)/max(real(trace(C))/Nr,eps); %#ok<AGROW>
        end
    end,end
end
if isempty(ratios),normalizedMedian=NaN;normalizedP90=NaN;
else
 normalizedMedian=median(ratios);
 normalizedP90=empiricalQuantile(ratios,.9);
end
cal=struct('developmentGroups',developmentGroups(:)', ...
    'generatorSha256',chasi_sha256_file(which('chasi_v3_calibrate_strength')), ...
    'gainErrorScale',max(median(gains),eps), ...
    'normalizedErrorMedian',normalizedMedian, ...
    'normalizedErrorP90',normalizedP90, ...
    'scope','development data only; never recompute on independent seed', ...
    'unit','v6 gainErrorScale is mean TX-column squared channel norm; legacy modes keep total TX norm');
end

function value=empiricalQuantile(x,q)
v=sort(double(x(:)));assert(~isempty(v) && q>=0 && q<=1);
position=1+(numel(v)-1)*q;lo=floor(position);hi=ceil(position);
value=v(lo)+(position-lo)*(v(hi)-v(lo));
end
