function [BIDS] = GTMSeg(BIDS)
% This function runs gtmseg from FreeSurfer on all subjects defined in
% the BIDS structure
% Author: Martin Norgaard

% FreeSurfer GTMSeg
subject_list = bids.query(BIDS, 'subjects');

if BIDS.config.env.nproc > 1
    parpool('local',BIDS.config.env.nproc)
    parfor idx = 1:numel(subject_list)
        if ~exist(fullfile(BIDS.pth,'derivatives/freesurfer',['sub-' num2str(subject_list{idx}),'/mri/gtmseg.mgz']),'file')
        setenv('SUBJECTS_DIR',fullfile(BIDS.pth,'derivatives/freesurfer'))
        unix(['gtmseg --s sub-' num2str(subject_list{idx}) ' --xcerseg'])
        end
    end
    delete(gcp('nocreate'));
else
    for idx = 1:numel(subject_list)
        if ~exist(fullfile(BIDS.pth,'derivatives/freesurfer',['sub-' num2str(subject_list{idx}),'/mri/gtmseg.mgz']),'file')
        setenv('SUBJECTS_DIR',fullfile(BIDS.pth,'derivatives/freesurfer'))
        unix(['gtmseg --s sub-' num2str(subject_list{idx}) ' --xcerseg'])
        end
    end
end


