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

mkdir(fullfile(BIDS.pth, 'derivatives',BIDS.config.env.derivatives_dir));
cd(fullfile(BIDS.pth, 'derivatives',BIDS.config.env.derivatives_dir));
for i = 1:numel(bids.query(BIDS, 'subjects'))
    cd(fullfile(BIDS.pth, 'derivatives',BIDS.config.env.derivatives_dir))
    mkdir(BIDS.subjects(i).name)
    cd(BIDS.subjects(i).name)
    mkdir(BIDS.subjects(i).session)
    cd(BIDS.subjects(i).session)
    mkdir('anat')
    mkdir('pet')
    end
end



