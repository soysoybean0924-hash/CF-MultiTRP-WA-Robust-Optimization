function out=chasi_scale_target_group_for_source(structural,targetPrecoder)
%CHASI_SCALE_TARGET_GROUP_FOR_SOURCE Exact one-group signal/power mapping.
% Target effective channel is H_t*x_t with per-group power weight nu.
% Source effective channel is H_s*X_s/sqrt(rb), power is sum |X_s|^2.
% H_s=sqrt(rb/nu)*H_t, X_s=sqrt(nu)*x_t makes both equal. For one source
% group its rate weight is 1; multiply reported SE by nu for target rate.
% This maps ALREADY noise-normalized arrays, not physical SRS parameters.
Htrue=structural.Htrue; Hdesign=structural.Hdesign;
K=size(Htrue,1); L=size(Htrue,2);
Nt=size(Htrue{1,1,1},2); d=size(targetPrecoder,4);
assert(isequal([size(targetPrecoder,1),size(targetPrecoder,2), ...
    size(targetPrecoder,3)],[L,K,Nt]));
assert(d>=1 && all(isfinite(targetPrecoder(:))));
nu=double(structural.resourceWeightNu);
rb=double(structural.rbCount);
if isfield(structural,'txRbCount'), rb=double(structural.txRbCount); end
assert(isfinite(nu) && nu>0 && isfinite(rb) && rb>0);
channelScale=sqrt(rb/nu); precoderScale=sqrt(nu);
S=size(structural.serviceTRP,2);
X=cell(K,1);
for k=1:K
    X{k,1}=complex(zeros(Nt*S,d));
    for l=1:L
        block=reshape(targetPrecoder(l,k,:,:),Nt,d);
        servingPosition=find(structural.serviceTRP(k,:)==l);
        if isempty(servingPosition)
            if norm(block,'fro')>1e-12
                error('chasi:NonservicePrecoder', ...
                    'Target precoder has nonzero data on a non-serving TRP.');
            end
        else
            rows=(servingPosition-1)*Nt+(1:Nt);
            X{k,1}(rows,:)=precoderScale*block;
        end
    end
end
for l=1:L
    for k=1:K
        Htrue{k,l,1}=channelScale*Htrue{k,l,1};
        Hdesign{k,l,1}=channelScale*Hdesign{k,l,1};
    end
end
out=struct('Htrue',{Htrue},'Hdesign',{Hdesign},'X',{X}, ...
    'serviceTRP',structural.serviceTRP,'ueIds',structural.ueIds, ...
    'rb',rb,'resourceWeightNu',nu,'sourceRateToTarget',nu, ...
    'channelScale',channelScale,'precoderScale',precoderScale);
end
