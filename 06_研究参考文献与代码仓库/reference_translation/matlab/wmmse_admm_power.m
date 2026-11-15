function out = wmmse_admm_power(a,B,Pmax,mode,draws,maxOuter,maxInner)
% Scalar power translation of WMMSE_ADMM_Func.py / WMMSE_PF_ADMM.py.
% a[L,K], B[L,L,K,K] are statistics from the original H/Hhat pipeline.
% draws[L,K,maxOuter] must be exported from the Python run for paired tests.
% Caps are diagnostic guards: reaching a cap is not convergence.
[L,K] = size(a);
Pmax=double(Pmax); % scipy.io.savemat preserves Python integer scalar type
assert(isequal(size(B),[L,L,K,K]));
assert(isequal(size(draws),[L,K,maxOuter]));
assert(strcmp(mode,'sumSE') || strcmp(mode,'PF'));
mu = 0.1*sqrt(Pmax/K)*ones(L,K);
rho = 0.001;
if strcmp(mode,'PF'), delta=0.001; else, delta=0.01; end
trace = nan(maxOuter+1,1);
primal = nan(maxOuter,maxInner); dual = primal;
innerCount = zeros(maxOuter,1);
trace(1) = utility(mu,a,B,mode);
previous = 0; current = trace(1); reason = 'outer_limit'; outer = 0;
while outer<maxOuter && (current-previous)^2>delta
    outer=outer+1; previous=current;
    numerator = sum(mu.*a,1).^2;
    denominator = ones(1,K);
    for k=1:K
        for i=1:K
            denominator(k)=denominator(k)+mu(:,i).'*B(:,:,k,i)*mu(:,i);
        end
    end
    v=sqrt(numerator)./denominator;
    e=1-numerator./denominator;
    if strcmp(mode,'PF'), w=-1./(e.*log(e)); else, w=1./e; end
    A=zeros(L,L,K); c=zeros(L,K);
    for k=1:K
        A(:,:,k)=(rho/2)*eye(L);
        for i=1:K
            A(:,:,k)=A(:,:,k)+w(i)*v(i)^2*B(:,:,i,k);
        end
        c(:,k)=w(k)*v(k)*a(:,k);
    end
    q=mu.*(1+draws(:,:,outer)); g=zeros(L,K);
    for j=1:maxInner
        oldq=q;
        for k=1:K
            mu(:,k)=A(:,:,k)\(c(:,k)+(rho/2)*(q(:,k)+g(:,k)));
        end
        if outer==1 && j==1
            debugFirstMu=mu; debugV=v; debugE=e; debugW=w; debugFirstQ=q;
        end
        q=mu-g;
        rownorm=sqrt(sum(q.^2,2));
        q=q.*min(1,sqrt(Pmax)./max(rownorm,eps));
        if outer==1 && j==1, debugFirstQProjected=q; end
        g=q-mu+g;
        primal(outer,j)=norm(mu-q,'fro');
        dual(outer,j)=rho*norm(q-oldq,'fro');
        innerCount(outer)=j;
        if j>5 && primal(outer,j)/norm(mu,'fro')<=0.001
            break
        end
    end
    current=utility(mu,a,B,mode);
    trace(outer+1)=current;
end
if (current-previous)^2<=delta, reason='original_outer_stop'; end
out.mu=mu.'; out.objective=trace(1:outer+1);
out.primal=primal(1:outer,:); out.dual=dual(1:outer,:);
out.innerCount=innerCount(1:outer); out.maxPower=max(sum(mu.^2,2));
out.stopReason=reason;
out.debugFirstMu=debugFirstMu; out.debugFirstQ=debugFirstQ;
out.debugFirstQProjected=debugFirstQProjected;
out.debugV=debugV; out.debugE=debugE; out.debugW=debugW;
end

function value=utility(mu,a,B,mode)
[~,K]=size(a);
numerator=sum(mu.*a,1).^2;
denominator=ones(1,K);
for k=1:K
    for i=1:K
        denominator(k)=denominator(k)+mu(:,i).'*B(:,:,k,i)*mu(:,i);
    end
end
se=log2(1+numerator./(denominator-numerator));
if strcmp(mode,'PF'), value=sum(log(se)); else, value=sum(se); end
end
