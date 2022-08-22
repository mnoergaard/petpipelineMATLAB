function BIDS = ConvertFS2BIDS(BIDS)
%
%
%

output_dir = fullfile(BIDS.pth,'derivatives',BIDS.config.env.derivatives_dir);
fs_dir = fullfile(BIDS.pth,'derivatives','freesurfer');
setenv('SUBJECTS_DIR',fullfile(BIDS.pth,'derivatives','freesurfer'));

if BIDS.config.env.nproc > 1
    parpool('local',BIDS.config.env.nproc)
    parfor idx = 1:numel(BIDS.subjects)
        subj = BIDS.subjects(idx).name;
        ses = BIDS.subjects(idx).session;
        aparc_file = fullfile(fs_dir, subj, 'mri/aparc+aseg.mgz');
        aparc_file_bids = fullfile(output_dir, subj, ses, ...
            'anat', [subj '_' ses '_space-T1w_desc-aparcaseg_dseg.nii.gz']);

        aseg_file = fullfile(fs_dir, subj, 'mri/aseg.mgz');
        aseg_file_bids = fullfile(output_dir, subj, ses, ...
            'anat', [subj '_' ses '_space-T1w_desc-aseg_dseg.nii.gz']);

        gtm_file = fullfile(fs_dir, subj, 'mri/gtmseg.mgz');
        gtm_file_bids = fullfile(output_dir, subj, ses, ...
            'anat', [subj '_' ses '_space-T1w_desc-gtm_dseg.nii.gz']);

        brainmask_file = fullfile(fs_dir, subj, 'mri/brainmask.mgz');
        brainmask_file_bids = fullfile(output_dir, subj, ses, ...
            'anat', [subj '_' ses '_space-T1w_desc-brain_mask.nii.gz']);

        unix(['mri_convert '  aparc_file ' ' aparc_file_bids]);
        unix(['mri_convert '  aseg_file ' ' aseg_file_bids]);
        unix(['mri_convert '  gtm_file ' ' gtm_file_bids]);
        unix(['mri_convert ' brainmask_file ' ' brainmask_file_bids]);

    end
    delete(gcp('nocreate'));
else
    for idx = 1:numel(BIDS.subjects)
        subj = BIDS.subjects(idx).name;
        ses = BIDS.subjects(idx).session;
        aparc_file = fullfile(fs_dir, subj, 'mri/aparc+aseg.mgz');
        aparc_file_bids = fullfile(output_dir, subj, ses, ...
            'anat', [subj '_' ses '_space-T1w_desc-aparcaseg_dseg.nii.gz']);

        aseg_file = fullfile(fs_dir, subj, 'mri/aseg.mgz');
        aseg_file_bids = fullfile(output_dir, subj, ses, ...
            'anat', [subj '_' ses '_space-T1w_desc-aseg_dseg.nii.gz']);

        gtm_file = fullfile(fs_dir, subj, 'mri/gtmseg.mgz');
        gtm_file_bids = fullfile(output_dir, subj, ses, ...
            'anat', [subj '_' ses '_space-T1w_desc-gtm_dseg.nii.gz']);

                brainmask_file = fullfile(fs_dir, subj, 'mri/brainmask.mgz');
        brainmask_file_bids = fullfile(output_dir, subj, ses, ...
            'anat', [subj '_' ses '_space-T1w_desc-brain_mask.nii.gz']);

        unix(['mri_convert '  aparc_file ' ' aparc_file_bids]);
        unix(['mri_convert '  aseg_file ' ' aseg_file_bids]);
        unix(['mri_convert '  gtm_file ' ' gtm_file_bids]);
        unix(['mri_convert ' brainmask_file ' ' brainmask_file_bids]);
    end
end