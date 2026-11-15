function summary=summarize_srs_overlap_paired(root)
%SUMMARIZE_SRS_OVERLAP_PAIRED Read-only same-input checkpoint summary.
base=fullfile(root,'07_实验结果与进展','local_runs', ...
    'srs_overlap_jakes_v3_main');
inputPath=fullfile(root,'03_研究场景构造','inputs', ...
    'srs_overlap_jakes_v3_main','input.mat');
loadedInput=load(inputPath,'data'); resources=loadedInput.data.resources;
names={'ideal_pf','nominal_pf','robust_pf'};
summary=struct('branch',names,'outerSteps',num2cell(zeros(1,3)), ...
    'complete',num2cell(false(1,3)),'stopReason',repmat({''},1,3), ...
    'inputSha256',repmat({''},1,3),'trueSumSE',num2cell(nan(1,3)), ...
    'trueR5',num2cell(nan(1,3)),'truePF',num2cell(nan(1,3)), ...
    'maxPower',num2cell(nan(1,3)),'totalPower',num2cell(nan(1,3)), ...
    'maxPowerViolation',num2cell(nan(1,3)), ...
    'innerSolved',num2cell(false(1,3)), ...
    'admmPrimal',num2cell(nan(1,3)),'admmDual',num2cell(nan(1,3)), ...
    'epsPrimal',num2cell(nan(1,3)),'epsDual',num2cell(nan(1,3)), ...
    'innerCount',num2cell(zeros(1,3)), ...
    'relativeObjectiveChange',num2cell(nan(1,3)), ...
    'precoderChange',num2cell(nan(1,3)));
for i=1:3
 path=fullfile(base,names{i},'latest.mat');
 if ~isfile(path), summary(i).stopReason='missing_checkpoint'; continue; end
 loaded=load(path,'state'); s=loaded.state;
 summary(i).outerSteps=s.iter; summary(i).complete=s.complete;
 summary(i).stopReason=s.stopReason;
 summary(i).inputSha256=s.inputSha256;
 if s.iter>0
  power=zeros(size(resources(1).H,2),1);
  for g=1:numel(resources)
   power=power+double(resources(g).resourceWeightNu)* ...
       reshape(sum(abs(s.precoder{g}).^2,[2 3 4]),[],1);
  end
  summary(i).totalPower=sum(power);
 end
 if s.iter>0
  e=s.history(end);
  summary(i).trueSumSE=e.trueSumSE;
  summary(i).trueR5=e.trueRate5;
  summary(i).truePF=e.truePF;
  summary(i).maxPower=e.maxPower;
  summary(i).maxPowerViolation=max(0,e.maxPower-s.nativeCfg.budget);
  summary(i).innerSolved=e.innerSolved;
  summary(i).admmPrimal=e.admmPrimal;
  summary(i).admmDual=e.admmDual;
  summary(i).epsPrimal=e.epsPrimal;
  summary(i).epsDual=e.epsDual;
  summary(i).innerCount=e.innerCount;
  summary(i).relativeObjectiveChange=e.relativeObjectiveChange;
  summary(i).precoderChange=e.precoderChange;
 end
end
ids={summary.inputSha256}; ids=ids(~cellfun(@isempty,ids));
assert(numel(unique(ids))<=1,'Branch inputs differ');
end
