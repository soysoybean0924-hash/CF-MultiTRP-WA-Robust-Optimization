function se=evaluate_original_se(a,B,mu,Pmax,prelog)
% Direct MATLAB translation of SpectralEfficiencyDownlink.py.
[L,K]=size(a); assert(isequal(size(mu),[K,L]));
mu=mu*sqrt(Pmax)/max(vecnorm(mu,2,1));
se=zeros(K,1);
for k=1:K
    numerator=(a(:,k).'*mu(k,:).')^2;
    denominator=1-numerator;
    for i=1:K
        denominator=denominator+mu(i,:)*B(:,:,k,i)*mu(i,:).';
    end
    se(k)=prelog*log2(1+numerator/denominator);
end
end
