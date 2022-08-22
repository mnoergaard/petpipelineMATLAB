function BIDS = PETVol2Vol(BIDS)
%
%
%

output_dir = fullfile(BIDS.pth,'derivatives',BIDS.config.env.derivatives_dir);
setenv('SUBJECTS_DIR',fullfile(BIDS.pth,'derivatives','freesurfer'));

if BIDS.config.env.nproc > 1
    parpool('local',BIDS.config.env.nproc)
    parfor idx = 1:numel(BIDS.subjects)
        subj = BIDS.subjects(idx).name;
        ses = BIDS.subjects(idx).session;

        switch BIDS.config.preproc.pvc.pvc
            case {'nopvc', 'agtm', 'gtm'}
                    input_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'input.nii.gz');
                output_file = fullfile(output_dir, subj, ses, ...
                    'pet',[subj '_' ses '_space-mni305_pvc-nopvc_desc-preproc_pet.nii.gz']);
            case 'mgx'
                input_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'mgx.subctxgm.nii.gz');
                output_file = fullfile(output_dir, subj, ses, ...
                    'pet',[subj '_' ses '_space-mni305_pvc-mgx_desc-preproc_pet.nii.gz']);
            case 'rbv'
                input_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'rbv.subctxgm.nii.gz');
                output_file = fullfile(output_dir, subj, ses, ...
                    'pet',[subj '_' ses '_space-mni305_pvc-rbv_desc-preproc_pet.nii.gz']);
        end

        lta_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'aux/bbpet2anat.lta');

        unix(['mri_vol2vol --mov '  input_file ...
            ' --reg ' lta_file ...
            ' --tal' ...
            ' --talres 2' ...
            ' --o ' output_file]);
    end
    delete(gcp('nocreate'));
else
    for idx = 1:numel(BIDS.subjects)
        subj = BIDS.subjects(idx).name;
        ses = BIDS.subjects(idx).session;

        switch BIDS.config.preproc.pvc.pvc
            case {'nopvc', 'agtm', 'gtm'}
                    input_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'input.nii.gz');
                output_file = fullfile(output_dir, subj, ses, ...
                    'pet',[subj '_' ses '_space-mni305_pvc-nopvc_desc-preproc_pet.nii.gz']);
            case 'mgx'
                input_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'mgx.subctxgm.nii.gz');
                output_file = fullfile(output_dir, subj, ses, ...
                    'pet',[subj '_' ses '_space-mni305_pvc-mgx_desc-preproc_pet.nii.gz']);
            case 'rbv'
                input_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'rbv.subctxgm.nii.gz');
                output_file = fullfile(output_dir, subj, ses, ...
                    'pet',[subj '_' ses '_space-mni305_pvc-rbv_desc-preproc_pet.nii.gz']);
        end

        lta_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'aux/bbpet2anat.lta');

        unix(['mri_vol2vol --mov '  input_file ...
            ' --reg ' lta_file ...
            ' --tal' ...
            ' --talres 2' ...
            ' --o ' output_file]);
    end
end