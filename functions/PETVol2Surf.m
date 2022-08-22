function BIDS = PETVol2Surf(BIDS)
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

        hemi = {'rh', 'lh'};
        for n = 1:numel(hemi)

            switch BIDS.config.preproc.pvc.pvc
                case {'nopvc', 'agtm', 'gtm'}
                    input_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'input.nii.gz');
                    output_file = fullfile(output_dir, subj, ses, ...
                        'pet',[subj '_' ses '_space-fsaverage_pvc-nopvc_hemi-' hemi{n} '_pet.nii.gz']);
                case 'mgx'
                    input_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'mgx.ctxgm.nii.gz');
                    output_file = fullfile(output_dir, subj, ses, ...
                        'pet',[subj '_' ses '_space-fsaverage_pvc-mgx_hemi-' hemi{n} '_pet.nii.gz']);
                case 'rbv'
                    input_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'rbv.ctxgm.nii.gz');
                    output_file = fullfile(output_dir, subj, ses, ...
                        'pet',[subj '_' ses '_space-fsaverage_pvc-rbv_hemi-' hemi{n} '_pet.nii.gz']);
            end

            lta_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'aux/bbpet2anat.lta');

            unix(['mri_vol2surf --mov '  input_file ...
                ' --reg ' lta_file ...
                ' --hemi ' hemi{n} ...
                ' --projfrac 0.5' ...
                ' --cortex' ...
                ' --trgsubject fsaverage' ...
                ' --o ' output_file]);
        end
    end
    delete(gcp('nocreate'));
else
    for idx = 1:numel(BIDS.subjects)
        subj = BIDS.subjects(idx).name;
        ses = BIDS.subjects(idx).session;

        hemi = {'rh', 'lh'};
        for n = 1:numel(hemi)

            switch BIDS.config.preproc.pvc.pvc
                case {'nopvc', 'agtm', 'gtm'}
                    input_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'input.nii.gz');
                    output_file = fullfile(output_dir, subj, ses, ...
                        'pet',[subj '_' ses '_space-fsaverage_pvc-nopvc_hemi-' hemi{n} '_pet.nii.gz']);
                case 'mgx'
                    input_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'mgx.ctxgm.nii.gz');
                    output_file = fullfile(output_dir, subj, ses, ...
                        'pet',[subj '_' ses '_space-fsaverage_pvc-mgx_hemi-' hemi{n} '_pet.nii.gz']);
                case 'rbv'
                    input_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'rbv.ctxgm.nii.gz');
                    output_file = fullfile(output_dir, subj, ses, ...
                        'pet',[subj '_' ses '_space-fsaverage_pvc-rbv_hemi-' hemi{n} '_pet.nii.gz']);
            end

            lta_file = fullfile(output_dir, subj, ses, ...
                'pet', BIDS.config.preproc.pvc.pvc, 'aux/bbpet2anat.lta');

            unix(['mri_vol2surf --mov '  input_file ...
                ' --reg ' lta_file ...
                ' --hemi ' hemi{n} ...
                ' --projfrac 0.5' ...
                ' --cortex' ...
                ' --trgsubject fsaverage' ...
                ' --o ' output_file]);
        end
    end
end