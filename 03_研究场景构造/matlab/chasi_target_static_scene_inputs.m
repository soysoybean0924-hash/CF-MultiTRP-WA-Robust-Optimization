function [resources,diagnostics]=chasi_target_static_scene_inputs(scene,model)
%CHASI_TARGET_STATIC_SCENE_INPUTS Isolated static estimator on 17 groups.
% Retains the target scene's H_true, schedule, RB weights and service mask;
% replaces ONLY estimated H and error statistics. No age model is applied.
resources=scene.resources;
G=numel(resources);
diagnostics=repmat(struct('group',0,'ueCount',0,'usedUserPilots',0,'nmse',NaN, ...
    'pilotNoiseVariance',NaN,'referenceGain',NaN),1,G);
for g=1:G
    groupModel=model;
    groupModel.seed=model.seed+g-1;
    estimate=chasi_measure_target_group_static(resources(g),groupModel);
    resources(g).H=single(estimate.Hhat);
    resources(g).C=single(permute(estimate.Ce,[3,4,1,2]));
    resources(g).C_stat=resources(g).C;
    diagnostics(g)=struct('group',g,'ueCount',numel(resources(g).ueIds), ...
        'usedUserPilots',estimate.usedUserPilots, ...
        'nmse',estimate.globalNmse, ...
        'pilotNoiseVariance',estimate.pilotNoiseVariance, ...
        'referenceGain',estimate.referenceGain);
end
end
