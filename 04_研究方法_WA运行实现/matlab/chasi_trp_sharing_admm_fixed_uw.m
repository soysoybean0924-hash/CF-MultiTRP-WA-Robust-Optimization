function out=chasi_trp_sharing_admm_fixed_uw(resources,parts,initial,budget,opts)
%CHASI_TRP_SHARING_ADMM_FIXED_UW Sharing ADMM on the same fixed-U/W V QP.
% f_l(V_l)=sum_g nu_g[D_gl||V_gl||^2-2Re<B_gl,V_gl>]+power-ball;
% g(sum_l A_l V_l)=sum_gk nu_g <S_gk,M_gk S_gk>,
% where (A_l V_l)_gk=H_gkl V_gl. Local x_l updates keep the across-group
% physical TRP power ball; the coordinator computes the aggregate prox.
% This MATLAB code simulates TRP messages serially. Parallel wall time is
% neither implemented nor inferred from this serial runtime.
arguments
 resources (1,:) struct
 parts (1,:) cell
 initial (1,:) cell
 budget (1,1) double {mustBePositive}
 opts.rho (1,1) double {mustBePositive}=1
 opts.maxIterations (1,1) double {mustBeInteger,mustBePositive}=400
 opts.absoluteTolerance (1,1) double {mustBePositive}=1e-7
 opts.relativeTolerance (1,1) double {mustBePositive}=1e-6
end
G=numel(resources); L=size(resources(1).H,2); assert(numel(parts)==G);
x=initial; rho=opts.rho; y=cell(G,L); u=cell(G,L);
for g=1:G
    [K,~,Nr,~]=size(resources(g).H); cols=K*size(x{g},4);
    for l=1:L
        y{g,l}=contribution(resources(g),x{g},l);
        u{g,l}=complex(zeros(Nr,cols,K));
    end
end
t0=tic; history=struct('iteration',{},'objective',{}, ...
    'primal',{},'dual',{},'epsPrimal',{},'epsDual',{}, ...
    'maxPowerViolation',{},'elapsedSeconds',{});
converged=false;
for iteration=1:opts.maxIterations
    for l=1:L
        Q=cell(1,G); rhs=cell(1,G); active=cell(1,G); nu=zeros(1,G);
        for g=1:G
            r=resources(g); [K,~,Nr,Nt]=size(r.H);
            streams=size(x{g},4); nu(g)=double(r.resourceWeightNu);
            active{g}=find(repmat(logical(r.servingMask(l,:)),1,streams));
            Q{g}=2*nu(g)*double(parts{g}.D(l))*eye(Nt);
            rhs{g}=2*nu(g)*block(parts{g}.B,l,K,Nt,streams);
            for k=1:K
                h=reshape(double(r.H(k,l,:,:)),Nr,Nt);
                target=y{g,l}(:,:,k)-u{g,l}(:,:,k);
                Q{g}=Q{g}+rho*(h'*h);
                rhs{g}=rhs{g}+rho*h'*target;
            end
            Q{g}=(Q{g}+Q{g}')/2;
        end
        candidate=localSolve(0,Q,rhs,active,nu);
        if weightedPower(candidate,nu)>budget
            low=0; high=1;
            while weightedPower(localSolve(high,Q,rhs,active,nu),nu)>budget
                high=2*high;
                assert(high<1e12,'Unable to bracket local ADMM power multiplier');
            end
            for j=1:35
                mid=(low+high)/2;
                if weightedPower(localSolve(mid,Q,rhs,active,nu),nu)>budget
                    low=mid;
                else
                    high=mid;
                end
            end
            candidate=localSolve(high,Q,rhs,active,nu);
        end
        for g=1:G
            K=size(resources(g).H,1); Nt=size(resources(g).H,4);
            streams=size(x{g},4);
            x{g}(l,:,:,:)=permute(reshape(candidate{g},Nt,K,streams), ...
                [4 2 1 3]);
        end
    end
    oldY=y; primalSquared=0; dualSquared=0; ySquared=0; aSquared=0;
    for g=1:G
        r=resources(g); [K,~,Nr,Nt]=size(r.H);
        cols=K*size(x{g},4); nu=double(r.resourceWeightNu);
        a=cell(1,L); total=complex(zeros(Nr,cols,K));
        for l=1:L
            a{l}=contribution(r,x{g},l);
            total=total+a{l}+u{g,l};
        end
        z=complex(zeros(size(total)));
        for k=1:K
            M=reshape(double(parts{g}.M(k,:,:)),Nr,Nr);
            z(:,:,k)=(2*nu*M+(rho/L)*eye(Nr))\ ...
                ((rho/L)*total(:,:,k));
        end
        shift=(z-total)/L;
        for l=1:L
            y{g,l}=a{l}+u{g,l}+shift;
            u{g,l}=u{g,l}+a{l}-y{g,l};
            primalSquared=primalSquared+sum(abs(a{l}-y{g,l}).^2,'all');
            aSquared=aSquared+sum(abs(a{l}).^2,'all');
            ySquared=ySquared+sum(abs(y{g,l}).^2,'all');
            deltaY=y{g,l}-oldY{g,l};
            adjoint=complex(zeros(Nt,cols));
            for k=1:K
                h=reshape(double(r.H(k,l,:,:)),Nr,Nt);
                adjoint=adjoint+h'*deltaY(:,:,k);
            end
            dualSquared=dualSquared+sum(abs(adjoint).^2,'all');
        end
    end
    primal=sqrt(primalSquared); dual=rho*sqrt(dualSquared);
    messageCount=0; variableCount=0;
    for g=1:G
        [K,~,Nr,Nt]=size(resources(g).H); streams=size(x{g},4);
        messageCount=messageCount+L*Nr*K*K*streams;
        variableCount=variableCount+L*Nt*K*streams;
    end
    epsP=sqrt(messageCount)*opts.absoluteTolerance+ ...
        opts.relativeTolerance*max(sqrt(aSquared),sqrt(ySquared));
    dualNorm=0;
    for g=1:G
        [K,~,Nr,Nt]=size(resources(g).H); streams=size(x{g},4);
        for l=1:L
            adjoint=complex(zeros(Nt,K*streams));
            for k=1:K
                h=reshape(double(resources(g).H(k,l,:,:)),Nr,Nt);
                adjoint=adjoint+h'*u{g,l}(:,:,k);
            end
            dualNorm=dualNorm+sum(abs(adjoint).^2,'all');
        end
    end
    epsD=sqrt(variableCount)*opts.absoluteTolerance+ ...
        opts.relativeTolerance*rho*sqrt(dualNorm);
    [objective,power]=fixedObjective(resources,parts,x);
    history(iteration)=struct('iteration',iteration, ...
        'objective',objective,'primal',primal,'dual',dual, ...
        'epsPrimal',epsP,'epsDual',epsD, ...
        'maxPowerViolation',max(max(power-budget,0)), ...
        'elapsedSeconds',toc(t0)); %#ok<AGROW>
    if mod(iteration,100)==0
        fprintf('TRP_SHARING_ADMM iter=%d objective=%.8g primal=%.3g/%.3g dual=%.3g/%.3g\n', ...
            iteration,objective,primal,epsP,dual,epsD);
    end
    if primal<=epsP && dual<=epsD
        converged=true; break
    end
end
out=struct('precoder',{x},'history',history, ...
    'method','TRP_local_sharing_ADMM_serial_simulation', ...
    'converged',converged,'parallelRuntimeMeasured',false, ...
    'coordinationComplexScalarsPerIteration',messageCount, ...
    'messageDescription','each TRP sends A_l V_l contribution; coordinator returns y_l and scaled dual', ...
    'elapsedSeconds',toc(t0),'fixedObjective',objective, ...
    'power',power,'maxPowerViolation',max(max(power-budget,0)));
end

function A=contribution(r,x,l)
[K,~,Nr,Nt]=size(r.H); streams=size(x,4);
X=block(x,l,K,Nt,streams); A=complex(zeros(Nr,K*streams,K));
for k=1:K
    h=reshape(double(r.H(k,l,:,:)),Nr,Nt);
    A(:,:,k)=h*X;
end
end

function X=block(v,l,K,Nt,streams)
X=reshape(permute(v(l,:,:,:),[3 2 4 1]),Nt,K*streams);
end

function result=localSolve(lambda,Q,rhs,active,nu)
result=cell(size(Q));
for g=1:numel(Q)
    result{g}=complex(zeros(size(rhs{g})));
    cols=active{g}; if isempty(cols), continue; end
    result{g}(:,cols)=(Q{g}+2*nu(g)*lambda*eye(size(Q{g})))\ ...
        rhs{g}(:,cols);
end
end

function power=weightedPower(x,nu)
power=0;
for g=1:numel(x), power=power+nu(g)*sum(abs(x{g}).^2,'all'); end
end

function [value,power]=fixedObjective(resources,parts,x)
G=numel(resources); L=size(resources(1).H,2);
power=zeros(L,1); value=0;
for g=1:G
    r=resources(g); [K,~,Nr,Nt]=size(r.H);
    streams=size(x{g},4); nu=double(r.resourceWeightNu); local=0;
    for l=1:L
        X=block(x{g},l,K,Nt,streams);
        B=block(parts{g}.B,l,K,Nt,streams);
        normSquared=sum(abs(X).^2,'all');
        power(l)=power(l)+nu*normSquared;
        local=local+double(parts{g}.D(l))*normSquared- ...
            2*real(sum(conj(B).*X,'all'));
    end
    for k=1:K
        S=complex(zeros(Nr,K*streams));
        for l=1:L
            h=reshape(double(r.H(k,l,:,:)),Nr,Nt);
            S=S+h*block(x{g},l,K,Nt,streams);
        end
        M=reshape(double(parts{g}.M(k,:,:)),Nr,Nr);
        local=local+real(sum(conj(S).*(M*S),'all'));
    end
    value=value+nu*local;
end
end
