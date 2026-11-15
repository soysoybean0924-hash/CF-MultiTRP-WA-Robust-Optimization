function out=chasi_trp_local_fixed_uw(resources,parts,initial,budget,maxSweeps)
%CHASI_TRP_LOCAL_FIXED_UW Serial TRP-block reference for fixed U/W.
% Each TRP solves its own convex V block across all resource groups under
% its physical total-power budget. The coordinator supplies aggregate
% received signals S_gk=sum_l H_gkl V_gl. This is block coordinate descent,
% NOT distributed ADMM or a measured parallel implementation.
G=numel(resources); L=size(resources(1).H,2);
assert(numel(parts)==G && numel(initial)==G && budget>0);
v=initial; t0=tic;
history=struct('sweep',{},'fixedObjective',{}, ...
    'maxPowerViolation',{},'relativePrecoderChange',{}, ...
    'elapsedSeconds',{});
for sweep=1:maxSweeps
    old=v;
    for l=1:L
        Q=cell(1,G); rhs=cell(1,G); active=cell(1,G); nu=zeros(1,G);
        for g=1:G
            r=resources(g); [K,~,Nr,Nt]=size(r.H);
            rank=size(v{g},4); cols=K*rank;
            nu(g)=double(r.resourceWeightNu);
            active{g}=find(repmat(logical(r.servingMask(l,:)),1,rank));
            Q{g}=double(parts{g}.D(l))*eye(Nt);
            rhs{g}=matrixBlock(parts{g}.B,l,K,Nt,rank);
            local=matrixBlock(v{g},l,K,Nt,rank);
            for k=1:K
                M=reshape(double(parts{g}.M(k,:,:)),Nr,Nr);
                h=reshape(double(r.H(k,l,:,:)),Nr,Nt);
                total=complex(zeros(Nr,cols));
                for j=1:L
                    hj=reshape(double(r.H(k,j,:,:)),Nr,Nt);
                    total=total+hj*matrixBlock(v{g},j,K,Nt,rank);
                end
                Q{g}=Q{g}+h'*M*h;
                rhs{g}=rhs{g}-h'*M*(total-h*local);
            end
            Q{g}=(Q{g}+Q{g}')/2;
        end
        candidate=solveAt(0,Q,rhs,active);
        if weightedPower(candidate,nu)>budget
            lo=0; hi=1;
            while weightedPower(solveAt(hi,Q,rhs,active),nu)>budget
                hi=2*hi;
                assert(hi<1e12,'Unable to bracket local power multiplier');
            end
            for step=1:35
                mid=(lo+hi)/2;
                if weightedPower(solveAt(mid,Q,rhs,active),nu)>budget
                    lo=mid;
                else
                    hi=mid;
                end
            end
            candidate=solveAt(hi,Q,rhs,active);
        end
        for g=1:G
            K=size(resources(g).H,1); Nt=size(resources(g).H,4);
            rank=size(v{g},4);
            v{g}(l,:,:,:)=permute(reshape(candidate{g},Nt,K,rank), ...
                [4 2 1 3]);
        end
    end
    [objective,power]=fixedObjective(resources,parts,v);
    delta=0; normOld=0;
    for g=1:G
        delta=delta+sum(abs(v{g}-old{g}).^2,'all');
        normOld=normOld+sum(abs(old{g}).^2,'all');
    end
    change=sqrt(delta)/max(1,sqrt(normOld));
    history(sweep)=struct('sweep',sweep,'fixedObjective',objective, ...
        'maxPowerViolation',max(max(power-budget,0)), ...
        'relativePrecoderChange',change,'elapsedSeconds',toc(t0)); %#ok<AGROW>
    if sweep>1
        rel=abs(history(sweep).fixedObjective- ...
            history(sweep-1).fixedObjective)/ ...
            max(1,abs(history(sweep-1).fixedObjective));
        if rel<1e-7 && change<1e-5, break; end
    end
end
out=struct('precoder',{v},'history',history, ...
    'method','serial_TRP_block_coordinate_fixed_UW_not_ADMM', ...
    'parallelRuntimeMeasured',false, ...
    'coordinationMessage','per resource/target UE aggregate Nr-by-(K*rank) complex S', ...
    'coordinationComplexScalarsPerSweep',sum(arrayfun(@(r) ...
        size(r.H,1)^2*size(r.H,3)*size(v{1},4),resources)), ...
    'localInformation','own H_gkl, WMMSE M_gk and B_gl, current aggregate S_gk, shared power budget', ...
    'elapsedSeconds',toc(t0),'fixedObjective',objective, ...
    'power',power,'maxPowerViolation',max(max(power-budget,0)));
end

function X=matrixBlock(x,l,K,Nt,rank)
X=reshape(permute(x(l,:,:,:),[3 2 4 1]),Nt,K*rank);
end

function result=solveAt(lambda,Q,rhs,active)
result=cell(size(Q));
for g=1:numel(Q)
    result{g}=complex(zeros(size(rhs{g})));
    cols=active{g};
    if isempty(cols), continue; end
    result{g}(:,cols)=(Q{g}+lambda*eye(size(Q{g})))\rhs{g}(:,cols);
end
end

function p=weightedPower(x,nu)
p=0;
for g=1:numel(x), p=p+nu(g)*sum(abs(x{g}).^2,'all'); end
end

function [objective,power]=fixedObjective(resources,parts,v)
G=numel(resources); L=size(resources(1).H,2);
power=zeros(L,1); objective=0;
for g=1:G
    r=resources(g); [K,~,Nr,Nt]=size(r.H);
    rank=size(v{g},4); nu=double(r.resourceWeightNu);
    local=0;
    for l=1:L
        X=matrixBlock(v{g},l,K,Nt,rank);
        B=matrixBlock(parts{g}.B,l,K,Nt,rank);
        energy=sum(abs(X).^2,'all');
        power(l)=power(l)+nu*energy;
        local=local+double(parts{g}.D(l))*energy- ...
            2*real(sum(conj(B).*X,'all'));
    end
    for k=1:K
        S=complex(zeros(Nr,K*rank));
        for l=1:L
            h=reshape(double(r.H(k,l,:,:)),Nr,Nt);
            S=S+h*matrixBlock(v{g},l,K,Nt,rank);
        end
        M=reshape(double(parts{g}.M(k,:,:)),Nr,Nr);
        local=local+real(sum(conj(S).*(M*S),'all'));
    end
    objective=objective+nu*local;
end
end
