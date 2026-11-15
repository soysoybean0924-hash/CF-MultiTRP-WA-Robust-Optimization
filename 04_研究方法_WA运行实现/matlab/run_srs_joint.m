function run_srs_joint(action,branch,root)
% Formal MATLAB entry; each call does bounded work, never launches background.
if nargin<3 || isempty(root), root=fileparts(fileparts(fileparts(mfilename('fullpath')))); end
if nargin<2 || isempty(branch), branch='all'; end
if nargin<1 || isempty(action), action='status'; end
addpath(root); startup();
configPath=fullfile(root,'03_研究场景构造','config','srs_joint_snapshot_v2.json');
cfg=jsondecode(fileread(configPath));
out=fullfile(root,'07_实验结果与进展','output_srs_joint_v2');
inputPath=fullfile(out,'input','srs_joint_snapshot_v2.mat');
names={'ideal_sumse','nominal_sumse','robust_sumse', ...
    'ideal_pf','nominal_pf','robust_pf'};
if strcmp(action,'input')
    assert(~isfile(inputPath),'Frozen input exists; choose a new output root.');
    prepare_srs_joint_input(root,configPath,inputPath);
    return
end
if strcmp(action,'export')
    export_srs_joint_binary(inputPath,fullfile(out,'frozen_binary'));
    return
end
if strcmp(action,'status')
    for i=1:numel(names)
        path=fullfile(out,'checkpoints',[names{i} '.mat']);
        if isfile(path)
            x=load(path,'state');
            fprintf('%s %d %d %s\n',names{i},x.state.iter, ...
                x.state.complete,x.state.stopReason);
        else, fprintf('%s absent\n',names{i}); end
    end
    return
end
if strcmp(action,'evaluate')
    assert(isfile(inputPath));
    dest=fullfile(out,'evaluation'); if ~isfolder(dest), mkdir(dest); end
    for i=1:numel(names)
        path=fullfile(out,'checkpoints',[names{i} '.mat']);
        if ~isfile(path), continue; end
        result=evaluate_srs_joint(inputPath,path); %#ok<NASGU>
        target=fullfile(dest,[names{i} '.mat']);
        assert(~isfile(target)); save(target,'result');
    end
    return
end
if strcmp(action,'plot')
    plot_srs_joint_zh(inputPath,fullfile(out,'checkpoints'),fullfile(out,'plots'));
    return
end
assert(strcmp(action,'solve') && isfile(inputPath));
if strcmp(branch,'all'), selected=names;
else, assert(ismember(branch,names)); selected={branch}; end
started=tic;
for i=1:numel(selected)
    if cfg.maxWallSeconds>0 && toc(started)>=cfg.maxWallSeconds
        fprintf('WALL_LIMIT_UNCONVERGED\n'); return
    end
    path=fullfile(out,'checkpoints',[selected{i} '.mat']);
    if isfile(path)
        x=load(path,'state'); if x.state.iter>=cfg.maxOuter
            fprintf('%s LIMIT_UNCONVERGED\n',selected{i}); continue
        end
    end
    run_srs_joint_batch(root,inputPath,selected{i},cfg.batchOuter,path);
end
end
