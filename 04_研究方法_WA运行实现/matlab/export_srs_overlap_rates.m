function result=export_srs_overlap_rates(root,csvPath)
%EXPORT_SRS_OVERLAP_RATES Per-UE static SE and distinct lower-tail mean.
base=fullfile(root,'07_实验结果与进展','local_runs','srs_overlap_jakes_v3_main');
names={'ideal_pf','nominal_pf','robust_pf'};
rate=zeros(210,3); tail=zeros(1,3); inputIds=cell(1,3);
for i=1:3
 loaded=load(fullfile(base,names{i},'latest.mat'),'state'); s=loaded.state;
 assert(s.iter>=1);
 rate(:,i)=double(s.history(end).truePerUserSE(:));
 inputIds{i}=s.inputSha256;
 sorted=sort(rate(:,i)); tail(i)=mean(sorted(1:ceil(.05*numel(sorted))));
end
assert(numel(unique(inputIds))==1);
writetable(table((1:210)',rate(:,1),rate(:,2),rate(:,3), ...
 'VariableNames',{'ueId','idealSE','nominalSE','totalCSE'}),csvPath);
result=struct('inputSha256',inputIds{1},'tail5Mean',tail, ...
 'tailDefinition','mean of lowest ceil(0.05*210) UE static SE values', ...
 'r5Definition','linearly interpolated 5th percentile');
end
