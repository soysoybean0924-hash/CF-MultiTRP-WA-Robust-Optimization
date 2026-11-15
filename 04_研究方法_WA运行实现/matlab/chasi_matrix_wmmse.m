function out = chasi_matrix_wmmse(resources, cfg, initialPrecoder)
%CHASI_MATRIX_WMMSE Complex coherent-JT WMMSE with central or TRP Jacobi ADMM.
% resources(g): H [K,L,Nr,Nt], C [K,L,Nr,Nr], servingMask [L,K],
%   rb_count, tx_rb_count. All inputs are supplied; no scene is generated.
% precoder{g} [L,K,Nt,numStreams], rxCombiner{g} [K,Nr,numStreams],
% mseWeight{g} and mseMatrix{g} [K,numStreams,numStreams].
% The local use of this function does not make its ADMM TRP distributed.
if nargin < 3, initialPrecoder = []; end
cfg = defaults(cfg);
wallClock=tic;
G = numel(resources); L = size(resources(1).H,2);
resources=prepareResources(resources,cfg);
if isempty(initialPrecoder)
    precoder = cell(1,G);
    for g=1:G
        resource=resources(g); [K,~,~,Nt]=size(resource.H);
        precoder{g}=complex(zeros(L,K,Nt,cfg.rank));
        for k=1:K
            for l=find(resource.servingMask(:,k)).'
                h=reshape(resource.H(k,l,:,:),size(resource.H,3),Nt);
                [~,~,vv]=svd(h,'econ');
                precoder{g}(l,k,:,:)=reshape(vv(:,1:cfg.rank),1,1,Nt,cfg.rank);
            end
        end
    end
    precoder=project(resources,precoder,cfg.budget);
    precoder=initializeVariation(resources,precoder,cfg);
else
    precoder=initialPrecoder;
    if numel(precoder)~=G, error('Resource count mismatch'); end
    for g=1:G
        resource=resources(g); expected=[L,size(resource.H,1),size(resource.H,4),cfg.rank];
        if ~isequal(size4(precoder{g}),expected), error('precoder shape mismatch'); end
        masked=precoder{g}.*reshape(~resource.servingMask,L,[],1,1);
        if sqrt(sum(abs(masked).^2,'all'))>1e-8
            error('precoder violates service mask');
        end
    end
    if max(chasi_power(resources,precoder))>cfg.budget+1e-7, error('precoder over budget'); end
end
hist=struct([]); rxCombiner=cell(1,G); mseWeight=cell(1,G); mseMatrix=cell(1,G); rate=cell(1,G);
lastObj=-inf; lastChange=inf; lastADMM=false; admmHist=[];
firstProbe=[];
for outer=0:cfg.outer_iterations
    for g=1:G
        [rxCombiner{g},mseWeight{g},mseMatrix{g},rate{g}]=state(resources(g).H,resources(g).C,precoder{g});
    end
    [obj,estimatedSum,perUE,userWeights,surrogate]=objectiveSnapshot(resources,rate,cfg);
    if outer==0, rel=NaN; else, rel=abs(obj-lastObj)/max(1,abs(lastObj)); end
    hist(outer+1).outer=outer; hist(outer+1).objective=obj;
    hist(outer+1).relative_objective_change=rel;
    hist(outer+1).max_power=max(chasi_power(resources,precoder));
    hist(outer+1).power_feasible= ...
        hist(outer+1).max_power<=cfg.budget+1e-6;
    hist(outer+1).service_feasible=serviceFeasible(resources,precoder);
    hist(outer+1).wall_seconds=toc(wallClock);
    hist(outer+1).stop_reason='continue';
    hist(outer+1).previous_precoder_change=lastChange;
    hist(outer+1).previous_admm_converged=lastADMM;
    hist(outer+1).estimated_sum_rate=estimatedSum;
    hist(outer+1).weighted_proxy=surrogate;
    if toc(wallClock)>cfg.max_wall_seconds
        reason='wall_budget_exceeded_at_outer_boundary';
        hist(outer+1).stop_reason=reason;break
    end
    if outer>0 && obj<lastObj-1e-7*max(1,abs(lastObj))
        reason='optimization_objective_decrease';
        hist(outer+1).stop_reason=reason;break
    end
    if outer>0 && rel<=cfg.outer_rel_tol && lastChange<=cfg.outer_precoder_tol && lastADMM
        reason='relative_objective_precoder_and_admm_tolerance';
        hist(outer+1).stop_reason=reason;break
    end
    if outer==cfg.outer_iterations
        if outer>0 && ~lastADMM, reason='max_outer_iterations_with_latest_inexact_admm';
        else, reason='max_outer_iterations'; end
        hist(outer+1).stop_reason=reason;break
    end
    parts=cell(1,G);
    for g=1:G
        % prepareResources has already added policyCov once to design C.
        parts{g}=quadratic(resources(g).H,resources(g).C,rxCombiner{g},mseWeight{g},userWeights{g});
    end
    if outer==0 && cfg.capture_first_step
        firstProbe=struct('H',{arrayfun(@(r)r.H,resources,'UniformOutput',false)}, ...
            'C',{arrayfun(@(r)r.C,resources,'UniformOutput',false)}, ...
            'A',{rxCombiner},'Omega',{mseWeight},'E',{mseMatrix}, ...
            'rate',{rate},'userWeights',{userWeights},'parts',{parts}, ...
            'objective',obj,'estimatedSum',estimatedSum,'perUE',perUE, ...
            'initialPrecoder',{precoder});
    end
    cfg.wallClockStart=wallClock;
    [candidatePrecoder,admmHist,lastADMM]=precoder_admm(resources,parts,precoder,cfg);
    if outer==0 && cfg.capture_first_step
        firstProbe.admmCandidate=candidatePrecoder;
        firstProbe.admmHistory=admmHist;
    end
    accepted=false; step=1; trialObj=-inf;
    for ls=1:12
        nextPrecoder=cell(1,G); trialRate=cell(1,G);
        for g=1:G
            nextPrecoder{g}=precoder{g}+step*(candidatePrecoder{g}-precoder{g});
            [~,~,~,trialRate{g}]=state(resources(g).H,resources(g).C,nextPrecoder{g});
        end
        trialObj=objectiveSnapshot(resources,trialRate,cfg);
        if trialObj>=obj-1e-7*max(1,abs(obj))
            accepted=true; break
        end
        step=step/2;
    end
    if accepted
        lastChange=cellnorm(cellminus(nextPrecoder,precoder))/max(1,cellnorm(precoder));
        precoder=nextPrecoder;
    else
        lastChange=0; step=0;
    end
    lastObj=obj;
    hist(outer+1).precoder_change=lastChange;
    hist(outer+1).admm_converged=lastADMM;
    hist(outer+1).admm_last=admmHist(end);
    hist(outer+1).admm_trace=admmHist;
    hist(outer+1).accepted_step=step;
    hist(outer+1).accepted_objective=trialObj;
    hist(outer+1).wall_seconds=toc(wallClock);
    if ~accepted
        reason='inexact_precoder_step_no_proxy_improvement';
        hist(outer+1).stop_reason=reason;break
    end
end
for g=1:G
    [rxCombiner{g},mseWeight{g},mseMatrix{g},rate{g}]=state(resources(g).H,resources(g).C,precoder{g});
end
[obj,estimatedSum,perUE,userWeights,surrogate]=objectiveSnapshot(resources,rate,cfg);
out.precoder=precoder; out.rxCombiner=rxCombiner; out.mseWeight=mseWeight;
out.mseMatrix=mseMatrix; out.rate=rate; out.trajectory=hist;
out.reason=reason; out.power=chasi_power(resources,precoder);
out.last_admm_history=admmHist;
out.optimization_objective=obj; out.estimated_sum_rate=estimatedSum;
out.solver_elapsed_seconds=toc(wallClock);
out.estimated_per_ue_rate=perUE; out.userWeight=userWeights;
out.weighted_proxy=surrogate; out.branch=cfg.branch;
out.converged=strcmp(reason,'relative_objective_precoder_and_admm_tolerance');
if ~strcmp(cfg.branch,'matrix_sum') && ...
        ~(isfield(cfg,'defer_true_evaluation') && cfg.defer_true_evaluation)
    out.true_per_ue_rate=trueRates(resources,precoder,cfg);
    out.true_rate5=linearPercentile(out.true_per_ue_rate,5);
    out.true_sum_rate=sum(out.true_per_ue_rate);
else
    out.true_per_ue_rate=[]; out.true_rate5=NaN; out.true_sum_rate=NaN;
end
out.max_power_violation=max(max(out.power-cfg.budget,0));
if cfg.capture_first_step, out.firstProbe=firstProbe; end
if cfg.probe_fixed_receiver_weight
    out.fixed_receiver_weight_probe=cell(1,G);
    for g=1:G
        coefficients=quadratic(resources(g).H,resources(g).C,rxCombiner{g},mseWeight{g},userWeights{g});
        qv=qapply(resources(g).H,coefficients,precoder{g},resources(g).servingMask);
        out.fixed_receiver_weight_probe{g}=struct('M',coefficients.M,...
            'B',coefficients.B,'D',coefficients.D,...
            'quadraticAppliedToPrecoder',qv,...
            'weighted_quadratic',resourceWeight(resources(g))*(real(sum(conj(precoder{g}).*qv,'all'))-...
              2*real(sum(conj(coefficients.B).*precoder{g},'all'))));
    end
end
end

function cfg=defaults(cfg)
defaultCfg=struct('rank',1,'budget',20,'outer_iterations',10,'outer_rel_tol',1e-4,...
    'outer_precoder_tol',5e-3,'admm_iterations',100,'cg_iterations',100,...
    'eps_abs',1e-4,'eps_rel',1e-3,'rho',.05,...
    'branch','matrix_sum','pf_floor',1e-2,'users',[],...
    'initializationMode','svd','initializationSeed',0,'seed',0,...
    'probe_fixed_receiver_weight',false,'max_wall_seconds',Inf);
defaultCfg.capture_first_step=false;
defaultCfg.coordinationMode='centralized';
f=fieldnames(defaultCfg);
for j=1:numel(f)
    if ~isfield(cfg,f{j}), cfg.(f{j})=defaultCfg.(f{j}); end
end
assert(ismember(char(cfg.coordinationMode), ...
    {'centralized','trp_jacobi_serial','trp_jacobi_parallel'}));
end

function precoder=initializeVariation(resources,precoder,cfg)
mode=char(cfg.initializationMode);
switch mode
    case 'svd'
        return
    case 'home_biased'
        for g=1:numel(resources)
            resource=resources(g); K=size(resource.H,1); L=size(resource.H,2);
            if ~isfield(resource,'home_trp') || numel(resource.home_trp)~=K
                error('Home-biased initialization requires home_trp');
            end
            home=double(resource.home_trp(:));
            if any(home<1 | home>L | home~=floor(home))
                error('Home-biased initialization requires one-based home_trp');
            end
            for k=1:K
                for l=find(resource.servingMask(:,k)).'
                    factor=.5;
                    if l==home(k), factor=1.5; end
                    precoder{g}(l,k,:,:)=factor*precoder{g}(l,k,:,:);
                end
            end
        end
    case 'jittered'
        seed=double(cfg.seed)+double(cfg.initializationSeed);
        if ~isscalar(seed) || ~isfinite(seed) || seed<0 || seed~=floor(seed)
            error('Jittered initialization requires a nonnegative integer seed');
        end
        stream=RandStream('mt19937ar','Seed',seed);
        for g=1:numel(resources)
            resource=resources(g); Nt=size(resource.H,4);
            noise=(randn(stream,size(precoder{g}))+1i*randn(stream,size(precoder{g})))/sqrt(2*Nt);
            mask=reshape(resource.servingMask,size(resource.servingMask,1),size(resource.servingMask,2),1,1);
            precoder{g}=precoder{g}+.15*noise.*mask;
        end
    otherwise
        error('Unsupported initialization mode: %s',mode);
end
precoder=project(resources,precoder,cfg.budget);
end

function resources=prepareResources(resources,cfg)
if isempty(resources), error('Native solver requires supplied resources'); end
if ~isfield(resources,'C_stat'), [resources.C_stat]=deal([]); end
if ~isfield(resources,'policyCov'), [resources.policyCov]=deal([]); end
if ~isfield(resources,'edge_weight'), [resources.edge_weight]=deal([]); end
L=size(resources(1).H,2); Nr=size(resources(1).H,3); Nt=size(resources(1).H,4);
for g=1:numel(resources)
    resource=resources(g); K=size(resource.H,1);
    if ~isequal(size4(resource.H),[K,L,Nr,Nt]) || ~isequal(size(resource.servingMask),[L,K])
        error('Resource %d has invalid H or serve axes',g);
    end
    if isempty(resource.C_stat)
        if isfield(resource,'C'), resource.C_stat=resource.C;
        else, error('Resource %d requires explicit error covariance',g); end
    end
    if ~isequal(size4(resource.C_stat),[K,L,Nr,Nr])
        error('Resource %d has invalid covariance axes',g);
    end
    nu=resourceWeight(resource);
    if ~isscalar(nu) || ~isfinite(nu) || nu<0 || nu>1
        error('Resource %d has invalid resourceWeightNu',g);
    end
    if ~strcmp(cfg.branch,'matrix_sum')
        if ~isfield(resource,'ueIds') || numel(resource.ueIds)~=K || ...
                any(resource.ueIds<1) || any(resource.ueIds~=floor(resource.ueIds))
            error('Resource %d requires one-based global ueIds',g);
        end
        if ismember(cfg.branch,{'ideal','ideal_pf'}) && ...
                (~isfield(resource,'H_true') || ...
                ~isequal(size4(resource.H_true),size4(resource.H)))
            error('Ideal branch requires explicit H_true',g);
        end
    end
    % Scene storage may be single precision. MATLAB's mixed single/double
    % matrix products otherwise round H*w before the double-precision MSE
    % solve, so identical explicit inputs disagree with the Python oracle.
    resource.H=double(resource.H);
    resource.C_stat=double(resource.C_stat);
    if ismember(cfg.branch,{'ideal','ideal_pf'})
        resource.H_true=double(resource.H_true);
    end
    switch cfg.branch
        case 'matrix_sum'
            resource.C=resource.C_stat;
        case 'ideal'
            resource.H=resource.H_true;
            resource.C=complex(zeros(K,L,Nr,Nr));
        case 'ideal_pf'
            resource.H=resource.H_true;
            resource.C=complex(zeros(K,L,Nr,Nr));
        case 'matrix_nominal'
            resource.C=complex(zeros(K,L,Nr,Nr));
        case {'matrix_statistical','matrix_edge_weighted'}
            resource.C=resource.C_stat;
        case 'matrix_pf_nominal'
            resource.C=complex(zeros(K,L,Nr,Nr));
        case 'matrix_pf_statistical'
            resource.C=resource.C_stat;
        case 'matrix_pf_uniform'
            resource.C=uniformCovariance(resource.C_stat,[]);
        case 'matrix_pf_uniform_relative'
            if ~isfield(resource,'prior_gain') || ~isequal(size(resource.prior_gain),[K,L])
                error('Uniform relative branch requires K-by-L prior_gain');
            end
            resource.C=uniformCovariance(resource.C_stat,resource.prior_gain);
        otherwise
            error('Unsupported native branch: %s',cfg.branch);
    end
    if ~isempty(resource.policyCov)
        assert(ismember(cfg.branch,{'matrix_statistical','matrix_pf_statistical'}));
        assert(isequal(size4(resource.policyCov),[K,L,Nr,Nr]));
        resource.policyCov=double(resource.policyCov);
        resource.C=resource.C+resource.policyCov;
    end
    if strcmp(cfg.branch,'matrix_edge_weighted')
        if ~isfield(resource,'home_trp') || numel(resource.home_trp)~=K || ...
                ~isfield(resource,'prior_gain') || ~isequal(size(resource.prior_gain),[K,L])
            error('Edge branch requires home_trp and K-by-L prior_gain');
        end
        home=double(resource.home_trp(:));
        if any(home<1 | home>L | home~=floor(home))
            error('Edge branch requires one-based home_trp');
        end
        gain=resource.prior_gain(sub2ind([K,L],(1:K).',home));
        if any(~isfinite(gain) | gain<0)
            error('Edge branch requires nonnegative finite home prior_gain');
        end
        potential=log2(1+cfg.budget*Nt*gain);
        potential=max(potential,.01);
        resource.edge_weight=min(8,max(.25,median(potential)./potential));
    end
    resources(g)=resource;
end
if ~strcmp(cfg.branch,'matrix_sum')
    allIDs=[];
    for g=1:numel(resources), allIDs=[allIDs;resources(g).ueIds(:)]; end
    if isempty(allIDs) || isempty(cfg.users) || cfg.users<max(allIDs)
        error('cfg.users must cover every one-based global ue_id');
    end
end
if cfg.rank<1 || cfg.rank>min(Nr,Nt) || cfg.pf_floor<=0 || cfg.budget<=0
    error('Invalid rank, PF floor, or budget');
end
end

function C=uniformCovariance(original,gain)
[K,L,Nr,~]=size(original); C=complex(zeros(size(original)));
active=reshape(any(any(abs(original)>0,4),3),K,L);
if ~any(active(:)), return; end
meanC=complex(zeros(Nr,Nr));
for k=1:K
    for l=1:L
        if active(k,l), meanC=meanC+reshape(original(k,l,:,:),Nr,Nr); end
    end
end
if isempty(gain)
    meanC=meanC/nnz(active);
else
    if any(~isfinite(gain(active))) || any(gain(active)<=0)
        error('Active relative-error links need positive finite prior_gain');
    end
    meanC=meanC/sum(gain(active));
end
meanC=(meanC+meanC')/2;
for k=1:K
    for l=1:L
        if active(k,l)
            factor=1; if ~isempty(gain), factor=gain(k,l); end
            C(k,l,:,:)=reshape(factor*meanC,1,1,Nr,Nr);
        end
    end
end
end

function [obj,sumRate,perUE,weights,surrogate]=objectiveSnapshot(resources,rate,cfg)
if strcmp(cfg.branch,'matrix_sum')
    obj=0; weights=cell(1,numel(resources));
    for g=1:numel(resources)
        obj=obj+resourceWeight(resources(g))*sum(rate{g});
        weights{g}=ones(numel(rate{g}),1);
    end
    sumRate=obj; perUE=[]; surrogate=obj; return
end
perUE=zeros(cfg.users,1); scheduled=false(cfg.users,1);
for g=1:numel(resources)
    ids=double(resources(g).ueIds(:));
    perUE(ids)=perUE(ids)+resourceWeight(resources(g))*rate{g}(:);
    scheduled(ids)=true;
end
sumRate=sum(perUE);
weights=cell(1,numel(resources));
switch cfg.branch
    case {'ideal_pf','matrix_pf_nominal','matrix_pf_statistical', ...
            'matrix_pf_uniform','matrix_pf_uniform_relative'}
        obj=sum(log1p(perUE(scheduled)/cfg.pf_floor));
        derivative=1./max(perUE+cfg.pf_floor,cfg.pf_floor);
        scale=median(derivative(scheduled));
        derivative=derivative/max(scale,1e-30);
        for g=1:numel(resources)
            weights{g}=derivative(double(resources(g).ueIds(:)));
        end
    case 'matrix_edge_weighted'
        obj=0;
        for g=1:numel(resources)
            weights{g}=resources(g).edge_weight(:);
            obj=obj+resourceWeight(resources(g))*sum(weights{g}.*rate{g}(:));
        end
    otherwise
        obj=sumRate;
        for g=1:numel(resources)
            weights{g}=ones(numel(rate{g}),1);
        end
end
surrogate=0;
for g=1:numel(resources)
    surrogate=surrogate+resourceWeight(resources(g))*sum(weights{g}.*rate{g}(:));
end
end

function perUE=trueRates(resources,precoder,cfg)
perUE=zeros(cfg.users,1);
for g=1:numel(resources)
    H=resources(g).H_true; x=precoder{g}; [K,L,Nr,Nt]=size(H); numStreams=size(x,4);
    for k=1:K
        allCov=eye(Nr); interf=eye(Nr);
        for i=1:K
            eff=complex(zeros(Nr,numStreams));
            for l=1:L
                eff=eff+reshape(H(k,l,:,:),Nr,Nt)*reshape(x(l,i,:,:),Nt,numStreams);
            end
            allCov=allCov+eff*eff';
            if i~=k, interf=interf+eff*eff'; end
        end
        spectralEfficiency=max(0,(real(log(det(allCov)))-real(log(det(interf))))/log(2));
        id=double(resources(g).ueIds(k));
        perUE(id)=perUE(id)+resourceWeight(resources(g))*spectralEfficiency;
    end
end
end

function value=linearPercentile(x,percentile)
x=sort(x(:)); n=numel(x);
if n==0, value=NaN; return; end
t=(n-1)*percentile/100; lo=floor(t)+1; hi=ceil(t)+1;
value=x(lo)+(t-floor(t))*(x(hi)-x(lo));
end

function [rxCombiner,mseWeight,mseMatrix,rate]=state(H,C,precoder)
[K,L,Nr,Nt]=size(H); numStreams=size(precoder,4);
rxCombiner=complex(zeros(K,Nr,numStreams)); mseWeight=complex(zeros(K,numStreams,numStreams));
mseMatrix=complex(zeros(K,numStreams,numStreams)); rate=zeros(K,1);
p=zeros(L,1);
for l=1:L, p(l)=sum(abs(precoder(l,:,:,:)).^2,'all'); end
for k=1:K
    rxCovariance=eye(Nr); wanted=complex(zeros(Nr,numStreams));
    for i=1:K
        effectiveChannel=complex(zeros(Nr,numStreams));
        for l=1:L
            h=reshape(H(k,l,:,:),Nr,Nt);
            precoderBlock=reshape(precoder(l,i,:,:),Nt,numStreams);
            effectiveChannel=effectiveChannel+h*precoderBlock;
        end
        rxCovariance=rxCovariance+effectiveChannel*effectiveChannel';
        if i==k, wanted=effectiveChannel; end
    end
    for l=1:L, rxCovariance=rxCovariance+p(l)*reshape(C(k,l,:,:),Nr,Nr); end
    combinerBlock=rxCovariance\wanted;
    mseBlock=eye(numStreams)-combinerBlock'*wanted-wanted'*combinerBlock+combinerBlock'*rxCovariance*combinerBlock;
    mseBlock=(mseBlock+mseBlock')/2;
    rxCombiner(k,:,:)=reshape(combinerBlock,1,Nr,numStreams);
    mseMatrix(k,:,:)=reshape(mseBlock,1,numStreams,numStreams);
    mseWeight(k,:,:)=reshape(mseBlock\eye(numStreams),1,numStreams,numStreams);
    rate(k)=-real(log(det(mseBlock)))/log(2);
end
end

function parts=quadratic(H,C,rxCombiner,mseWeight,userWeights)
[K,L,Nr,Nt]=size(H); numStreams=size(rxCombiner,3);
M=complex(zeros(K,Nr,Nr)); B=complex(zeros(L,K,Nt,numStreams)); D=zeros(L,1);
for k=1:K
    combinerBlock=reshape(rxCombiner(k,:,:),Nr,numStreams);
    weightedMse=reshape(mseWeight(k,:,:),numStreams,numStreams)*userWeights(k);
    quadraticBlock=combinerBlock*weightedMse*combinerBlock'; M(k,:,:)=reshape(quadraticBlock,1,Nr,Nr);
    for l=1:L
        h=reshape(H(k,l,:,:),Nr,Nt);
        B(l,k,:,:)=reshape(h'*combinerBlock*weightedMse,1,1,Nt,numStreams);
        c=reshape(C(k,l,:,:),Nr,Nr);
        D(l)=D(l)+real(trace(quadraticBlock*c));
    end
end
parts=struct('M',M,'B',B,'D',D);
end

function y=qapply(H,parts,x,mask)
[K,L,Nr,Nt]=size(H); numStreams=size(x,4); y=complex(zeros(size(x)));
active=cell(L,1); columns=cell(L,1); activeX=cell(L,1); activeY=cell(L,1);
for l=1:L
    active{l}=find(mask(l,:));
    columns{l}=reshape(active{l}(:)+(0:numStreams-1)*K,1,[]);
    activeX{l}=reshape(permute(x(l,active{l},:,:),[3 2 4 1]),Nt,[]);
    activeY{l}=complex(zeros(size(activeX{l})));
end
for k=1:K
    m=reshape(parts.M(k,:,:),Nr,Nr);
    % Stack every user's streams into one GEMM. This is the same quadratic
    % operator as the per-user loop, with the reduction order changed.
    s=complex(zeros(Nr,K*numStreams));
    for l=1:L
        if isempty(active{l}), continue; end
        h=reshape(H(k,l,:,:),Nr,Nt);
        cols=columns{l};
        s(:,cols)=s(:,cols)+h*activeX{l};
    end
    ms=m*s;
    for l=1:L
        if isempty(active{l}), continue; end
        h=reshape(H(k,l,:,:),Nr,Nt);
        activeY{l}=activeY{l}+h'*ms(:,columns{l});
    end
end
for l=1:L
    if isempty(active{l}), continue; end
    y(l,active{l},:,:)=permute(reshape(activeY{l},Nt,numel(active{l}),numStreams),[4 2 1 3]);
end
for l=1:L, y(l,:,:,:)=y(l,:,:,:)+parts.D(l)*x(l,:,:,:); end
y=y.*reshape(mask,L,K,1,1);
end

function [z,history,converged]=precoder_admm(resources,parts,initial,cfg)
G=numel(resources); z=initial; x=initial; dual=cellfun(@(precoderBlock) zeros(size(precoderBlock),'like',precoderBlock),z,'UniformOutput',false);
n=sum(cellfun(@numel,z)); rho=cfg.rho; converged=false; history=struct([]);
if strcmp(cfg.coordinationMode,'trp_jacobi_parallel')
    assert(license('test','Distrib_Computing_Toolbox') && ...
        ~isempty(gcp('nocreate')), ...
        'Parallel Computing Toolbox and an existing pool are required');
end
for t=1:cfg.admm_iterations
    new=cell(1,G); cgmax=0;
    exchangeBytes=0;coordinationRounds=0;
    if ~strcmp(cfg.coordinationMode,'centralized')
        coordinationRounds=1;
    end
    for g=1:G
        resource=resources(g); mask=reshape(resource.servingMask,size(resource.servingMask,1),size(resource.servingMask,2),1,1);
        nu=resourceWeight(resource); coefficients=parts{g};
        if strcmp(cfg.coordinationMode,'centralized')
            apply=@(trialPrecoder) nu*qapply(resource.H,coefficients,trialPrecoder,resource.servingMask)+(rho/2)*trialPrecoder;
            rhs=(nu*coefficients.B+(rho/2)*(z{g}-dual{g})).*mask;
            new{g}=cg(apply,rhs,x{g},cfg.cg_iterations,1e-5).*mask;
            residual=apply(new{g})-rhs;
            cgmax=max(cgmax,sqrt(sum(abs(residual).^2,'all'))/max(sqrt(sum(abs(rhs).^2,'all')),eps));
        else
            [new{g},relative,bytes]=trpJacobiUpdate( ...
                resource,coefficients,x{g},z{g},dual{g}, ...
                nu,rho,cfg.cg_iterations,cfg.coordinationMode);
            cgmax=max(cgmax,relative);
            exchangeBytes=exchangeBytes+bytes;
        end
    end
    oldz=z; z=project(resources,cellplus(new,dual),cfg.budget);
    dual=cellplus(dual,cellminus(new,z)); x=new;
    primal=cellnorm(cellminus(new,z)); dualres=rho*cellnorm(cellminus(z,oldz));
    ep=sqrt(n)*cfg.eps_abs+cfg.eps_rel*max(cellnorm(new),cellnorm(z));
    ed=sqrt(n)*cfg.eps_abs+cfg.eps_rel*rho*cellnorm(dual);
    history(t).inner=t;history(t).rho=rho;history(t).primal=primal;
    history(t).dual=dualres;history(t).eps_primal=ep;
    history(t).eps_dual=ed;history(t).cg_relative=cgmax;
    history(t).power=chasi_power(resources,z);
    history(t).coordinationMode=cfg.coordinationMode;
    history(t).coordinationRounds=coordinationRounds;
    history(t).estimatedExchangeBytes=exchangeBytes;
    history(t).exchangeScope= ...
        'logical complex-double receive aggregates to coordinator and back; not measured network traffic';
    history(t).synchronizationRule= ...
        'all TRP blocks use preceding z; consensus projection and dual update follow the complete round';
    history(t).max_power_violation=max(max(history(t).power-cfg.budget,0));
    history(t).wallSeconds=toc(cfg.wallClockStart);
    history(t).wallBudgetExceeded= ...
        history(t).wallSeconds>cfg.max_wall_seconds;
    if cfg.probe_fixed_receiver_weight
        fixedObjective=0;
        for g=1:G
            qz=qapply(resources(g).H,parts{g},z{g},resources(g).servingMask);
            fixedObjective=fixedObjective+resourceWeight(resources(g))*(...
                real(sum(conj(z{g}).*qz,'all'))-...
                2*real(sum(conj(parts{g}.B).*z{g},'all')));
        end
        history(t).fixed_receiver_weight_precoder_objective=fixedObjective;
    end
    if primal<=ep && dualres<=ed && cgmax<=max(1e-3,10*cfg.eps_rel)
        converged=true; break
    end
    if history(t).wallBudgetExceeded,break;end
    if mod(t,10)==0 && primal>10*dualres && rho<10
        oldrho=rho;rho=2*rho;dual=cellfun(@(dualBlock) dualBlock*oldrho/rho,dual,'UniformOutput',false);
    elseif mod(t,10)==0 && dualres>10*primal && rho>1e-5
        oldrho=rho;rho=rho/2;dual=cellfun(@(dualBlock) dualBlock*oldrho/rho,dual,'UniformOutput',false);
    end
end
end

function [updated,relative,bytes]=trpJacobiUpdate(resource,parts,x,z,dual,nu,rho,cgIterations,mode)
% One synchronous block-Jacobi round for fixed receive filters/weights.
% Each TRP solves its local diagonal Hessian against the previous round's
% aggregate. The coordinator forms that aggregate before all updates.
L=size(resource.H,2);
aggregateGradient=nu*qapply(resource.H,parts,z,resource.servingMask);
updated=complex(zeros(size(z)));residual=zeros(L,1);
if strcmp(mode,'trp_jacobi_parallel')
    blocks=cell(1,L);errors=zeros(1,L);
    parfor l=1:L
        [blocks{l},errors(l)]=trpLocalBlock( ...
            resource,parts,x,z,dual,aggregateGradient,nu,rho,cgIterations,l);
    end
    for l=1:L,updated(l,:,:,:)=blocks{l};end
    residual=errors(:);
else
    for l=1:L
        [updated(l,:,:,:),residual(l)]=trpLocalBlock( ...
            resource,parts,x,z,dual,aggregateGradient,nu,rho,cgIterations,l);
    end
end
relative=max(residual);
% Logical payload: each TRP's complex-double z block to the central
% aggregate, then the matching aggregate-gradient block back. This omits
% headers, transport, CSI distribution, and consensus projection traffic.
% Serial mode computes the same update in one process.
bytes=16*(numel(z)+numel(aggregateGradient));
end

function [block,relative]=trpLocalBlock(resource,parts,x,z,dual,aggregateGradient,nu,rho,cgIterations,l)
shape=size(z);L=shape(1);mask=resource.servingMask;
localZero=complex(zeros(shape));localZero(l,:,:,:)=z(l,:,:,:);
own=nu*qapply(resource.H,parts,localZero,mask);
cross=aggregateGradient(l,:,:,:)-own(l,:,:,:);
rhs=(nu*parts.B(l,:,:,:)-cross+ ...
    (rho/2)*(z(l,:,:,:)-dual(l,:,:,:))).* ...
    reshape(mask(l,:),1,[],1,1);
apply=@(trial)trpLocalApply(resource,parts,trial,l,nu,rho,shape);
block=cg(apply,rhs,x(l,:,:,:),cgIterations,1e-5).* ...
    reshape(mask(l,:),1,[],1,1);
error=apply(block)-rhs;
relative=sqrt(sum(abs(error).^2,'all'))/ ...
    max(sqrt(sum(abs(rhs).^2,'all')),eps);
end

function value=trpLocalApply(resource,parts,block,l,nu,rho,shape)
trial=complex(zeros(shape));trial(l,:,:,:)=block;
applied=qapply(resource.H,parts,trial,resource.servingMask);
value=nu*applied(l,:,:,:)+(rho/2)*block;
end

function x=cg(apply,b,x,maxiter,tol)
cgResidual=b-apply(x); cgDirection=cgResidual;
rr=real(sum(conj(cgResidual).*cgResidual,'all'));
% Match the source Python CG's relative tolerance, rtol*norm(b).
target=tol^2*real(sum(conj(b).*b,'all'));
for j=1:maxiter
    if rr<=target, break; end
    ap=apply(cgDirection);
    cgStepSize=rr/max(real(sum(conj(cgDirection).*ap,'all')),1e-30);
    x=x+cgStepSize*cgDirection; cgResidual=cgResidual-cgStepSize*ap;
    rrnew=real(sum(conj(cgResidual).*cgResidual,'all'));
    cgDirection=cgResidual+(rrnew/max(rr,1e-30))*cgDirection; rr=rrnew;
end
end

function out=project(resources,x,budget)
out=x; G=numel(x); L=size(x{1},1);
for g=1:G
    resource=resources(g); out{g}=out{g}.*reshape(resource.servingMask,L,[],1,1);
end
p=chasi_power(resources,out);
for l=1:L
    if p(l)<=budget, continue; end
    lo=0; hi=1;
    while projectedpower(resources,out,l,hi)>budget, hi=2*hi; end
    for j=1:35
        lam=(lo+hi)/2;
        if projectedpower(resources,out,l,lam)>budget, lo=lam;
        else, hi=lam; end
    end
    for g=1:G, out{g}(l,:,:,:)=out{g}(l,:,:,:)/(1+hi*resourceWeight(resources(g))); end
end
end

function p=projectedpower(resources,x,l,lambda)
p=0;
for g=1:numel(x)
    nu=resourceWeight(resources(g));
    p=p+nu*sum(abs(x{g}(l,:,:,:)).^2,'all')/(1+lambda*nu)^2;
end
end

function p=chasi_power(resources,precoder)
L=size(precoder{1},1); p=zeros(L,1);
for g=1:numel(precoder)
    nu=resourceWeight(resources(g));
    for l=1:L, p(l)=p(l)+nu*sum(abs(precoder{g}(l,:,:,:)).^2,'all'); end
end
end

function feasible=serviceFeasible(resources,precoder)
feasible=true;
for g=1:numel(resources)
    mask=resources(g).servingMask;
    outside=precoder{g}.*reshape(~mask,size(mask,1),size(mask,2),1,1);
    feasible=feasible && all(isfinite(precoder{g}),'all') && ...
        sum(abs(outside).^2,'all')<=1e-12;
end
end

function nu=resourceWeight(resource)
if isfield(resource,'resourceWeightNu') && ~isempty(resource.resourceWeightNu)
    nu=double(resource.resourceWeightNu);
elseif isfield(resource,'tx_rb_count')
    nu=double(resource.tx_rb_count)/273;
else
    nu=double(resource.rb_count)/273;
end
end

function s=size4(array)
s=[size(array,1),size(array,2),size(array,3),size(array,4)];
end

function n=cellnorm(c)
n=sqrt(sum(cellfun(@(x) real(sum(conj(x).*x,'all')),c)));
end

function c=cellplus(a,b)
c=cellfun(@plus,a,b,'UniformOutput',false);
end

function c=cellminus(a,b)
c=cellfun(@minus,a,b,'UniformOutput',false);
end

