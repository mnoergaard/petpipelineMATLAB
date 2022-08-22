function [] = create_dirs_derivative(BIDS, config)
% This function creates the derivatives directories for all subjects,
% sessions and modalities. A pipeline name is also requested to that will
% be located inside the derivatives directory.
% Author: Martin Norgaard, 01/26/2022

config_num = regexp(config,'\d*','Match');
if isempty(config_num)
    pipeline_name = BIDS.config.env.derivatives_dir;
else
    pipeline_name = [BIDS.config.env.derivatives_dir config_num{1}];
    BIDS.config.env.derivatives_dir = [BIDS.config.env.derivatives_dir config_num{1}];
end

mkdir(fullfile(BIDS.pth, 'derivatives'));
cd(fullfile(BIDS.pth, 'derivatives'));
for i = 1:numel(bids.query(BIDS, 'subjects'))
    subject_list = bids.query(BIDS, 'subjects');
    mkdir(fullfile(pipeline_name, ['sub-' num2str(subject_list{i})]));
    for j = 1:numel(bids.query(BIDS, 'sessions'))
        session_list = bids.query(BIDS, 'sessions');
        mkdir(fullfile(pipeline_name, ['sub-' num2str(subject_list{i})], ...
            ['ses-' session_list{j}]));
        for k = 1:numel(bids.query(BIDS, 'modalities'))
            modalities_list = bids.query(BIDS, 'modalities');
            mkdir(fullfile(pipeline_name, ['sub-' num2str(subject_list{i})], ...
                ['ses-' session_list{j}], ...
                modalities_list{k}));
        end
    end
end



