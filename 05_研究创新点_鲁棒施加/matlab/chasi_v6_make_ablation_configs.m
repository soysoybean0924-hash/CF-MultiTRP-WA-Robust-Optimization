function dir=chasi_v6_make_ablation_configs(root,inputPath,mode,strengthScale,developmentInput,proxyCalibrationPath)
% Lock development-group calibration once; every ablation is a separate
% six-branch campaign using the same input and four fixed baselines.
allowed={'source_all','source_service','source_beam','source_power', ...
    'source_residual','source_pollution','source_noise','source_age', ...
    'source_uniform_strength','source_position_strength', ...
    'source_csi_noise','source_power_constrained'};
assert(ismember(mode,allowed));startup;
if nargin<4,strengthScale=1;end
if nargin<5 || isempty(developmentInput),developmentInput=inputPath;end
if nargin<6,proxyCalibrationPath='';end
assert(isempty(proxyCalibrationPath) || ...
    strcmp(mode,'source_power_constrained'));
calibrationSha='';calibrationTag='uncalibrated';
if ~isempty(proxyCalibrationPath)
 assert(isfile(proxyCalibrationPath));
 calibrationSha=chasi_sha256_file(proxyCalibrationPath);
 calibrationTag=calibrationSha(1:16);
end
assert(isscalar(strengthScale) && isfinite(strengthScale) && ...
    strengthScale>=0 && strengthScale<=1 && ...
    abs(strengthScale*1000-round(strengthScale*1000))<1e-9);
inputSha=chasi_sha256_file(inputPath);
developmentSha=chasi_sha256_file(developmentInput);
assert(strcmp(inputSha,strtrim(fileread([inputPath '.sha256']))));
assert(strcmp(developmentSha, ...
    strtrim(fileread([developmentInput '.sha256']))));
base=fullfile(root,'05_研究创新点_鲁棒施加','config', ...
    ['v6_development_locked_' developmentSha(1:16)]);
if ~isfile(fullfile(base,'source_all.json'))
 assert(~isfolder(base) || isempty(builtin('dir',fullfile(base,'*.json'))), ...
     'Incomplete locked calibration; inspect before regenerating');
 chasi_v3_make_policy_configs(root,base,developmentInput,1:8);
end
locked=jsondecode(fileread(fullfile(base,'source_all.json')));
assert(isfield(locked.policy.calibration,'generatorSha256') && ...
    strcmp(locked.policy.calibration.generatorSha256, ...
        chasi_sha256_file(which('chasi_v3_calibrate_strength'))), ...
    'Development calibration implementation changed; use a new v6 input identity');
scaleTag='';
if strengthScale~=1,scaleTag=sprintf('_strength%04d',round(strengthScale*1000));end
dir=fullfile(root,'05_研究创新点_鲁棒施加','config', ...
    ['v6_' mode scaleTag '_' inputSha(1:16) '_dev' developmentSha(1:16)]);
if strcmp(mode,'source_power_constrained')
 dir=[dir '_proxy' calibrationTag];
end
if ~isfolder(dir),mkdir(dir);end
names={'nominal','physical','uniform'};
for j=1:numel(names)
 dest=fullfile(dir,[names{j} '.json']);
 if ~isfile(dest),copyfile(fullfile(base,[names{j} '.json']),dest);end
end
dest=fullfile(dir,'source_all.json');
if ~isfile(dest)
 if strcmp(mode,'source_uniform_strength') || ...
        strcmp(mode,'source_position_strength')
  s=jsondecode(fileread(fullfile(base,'source_all.json')));
  if strcmp(mode,'source_uniform_strength')
   s.policy.strengthMode='uniform';
  else
   s.policy.strengthMode='position';
   % Prespecified, bounded position contrast; no experimental tuning.
   s.policy.strength(:,1)=.05;
   s.policy.strength(:,2)=.10;
   s.policy.strength(:,3)=.15;
  end
  writeV6(dest,s,strengthScale,developmentSha);
 else
  if strcmp(mode,'source_power_constrained')
   s=jsondecode(fileread(fullfile(base,'source_all.json')));
   s.experienceConstraint=struct('enabled',true, ...
       'modeId','source_fixed_direction_power_proxy_v1', ...
       'decisionInputSha256',inputSha, ...
       'decisionSceneSeed',loadSceneSeed(inputPath), ...
       'maximumMeanLoss',.10, ...
       'proxyModelId','causal_Hhat_physicalC_Gaussian_logdet_v1', ...
       'proxyCalibrationPath',proxyCalibrationPath, ...
       'proxyCalibrationSha256',calibrationSha);
   writeV6(dest,s,strengthScale,developmentSha);
  elseif strcmp(mode,'source_csi_noise')
   s=jsondecode(fileread(fullfile(base,'source_all.json')));
   s.policy.enabled=false(4,3);
   writeV6(dest,s,strengthScale,developmentSha);
  else
   s=jsondecode(fileread(fullfile(base,[mode '.json'])));
   writeV6(dest,s,strengthScale,developmentSha);
  end
 end
end
end

function seed=loadSceneSeed(path)
x=load(path,'data');seed=x.data.sceneSeed;
end

function writeV6(path,s,strengthScale,developmentSha)
% The v5 source statistics merged SRS and CSI-RS noise. V6 retains the
% independently stored fifth covariance and gives it its own beam switch.
s.policy.enabled=[logical(s.policy.enabled);false(1,3)];
s.policy.strength=[s.policy.strength;s.policy.strength(3,:)];
s.policy.strength=s.policy.strength*strengthScale;
s.developmentInputSha256=developmentSha;
if strcmp(s.mode,'source')
 if ~any(s.policy.enabled(:))
  s.policy.enabled(5,2)=true; % CSI-RS-only ablation.
 else
  s.policy.enabled(5,2)=s.policy.enabled(3,2);
 end
end
fid=fopen(path,'w');assert(fid>0);
fwrite(fid,jsonencode(s,PrettyPrint=true),'char');fclose(fid);
end
