function BIDS = CoReg(BIDS)
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
        mc_name = fullfile(output_dir, subj, ses, ...
            'pet', [subj '_' ses '_desc-mc_pet.nii.gz']);

        output_mean = fullfile(output_dir, subj, ses, ...
            'pet', [subj '_' ses '_desc-mc_mean.nii.gz']);

        unix(['mri_concat '  mc_name ...
            ' --mean' ...
            ' --o ' output_mean]);

        mc_mean = struct;
        mc_mean.Description = 'Mean PET image averaged across all time frames';
        mc_mean.Sources = erase(mc_name,pwd);
        mc_mean.SoftwareName = 'FreeSurfer-mri_concat';
        [~, mc_mean.SoftwareVersion] = unix('mri_concat --version');
        mc_mean.CommandLine = ['mri_concat '  erase(mc_name,pwd) ...
            ' --mean' ...
            ' --o ' erase(output_mean,pwd)];

        bids.util.jsonwrite(fullfile(output_dir, subj, ses, ...
            'pet', [subj '_' ses '_desc-mc_mean.json']),mc_mean);

        output_lta = fullfile(output_dir, subj, ses, ...
            'pet', [subj '_' ses '_from-pet_to-T1w_reg.lta']);

        unix(['mri_coreg --s '  subj ...
            ' --mov ' output_mean ...
            ' --reg ' output_lta]);

        brainmask_file_bids = fullfile(output_dir, subj, ses, ...
            'anat', [subj '_' ses '_space-T1w_desc-brain_mask.nii.gz']);


        % Write metadata to json
        coreg = struct;
        coreg.Description = 'Registration file to go from PET space to ANAT space';
        coreg.Sources = {erase(output_mean,pwd), erase(brainmask_file_bids,pwd)};
        coreg.ReferenceImage = fullfile(erase(output_dir,pwd), subj, ses, ...
            'anat', [subj '_' ses '_space-T1w_desc-brain_mask.nii.gz']);
        coreg.dof = BIDS.config.preproc.coreg.dof;
        coreg.QC = {};
        coreg.SoftwareName = 'FreeSurfer-mri_coreg';
        [~, coreg.SoftwareVersion] = unix('mri_coreg --version');
        coreg.CommandLine = ['mri_coreg --s '  subj ...
            ' --mov ' erase(output_mean,pwd) ...
            ' --reg ' erase(output_lta, pwd)];
        
        bids.util.jsonwrite(fullfile(output_dir, subj, ses, ...
            'pet', [subj '_' ses '_from-pet_to-T1w_reg.json']),coreg);

    end
    delete(gcp('nocreate'));
else
    for idx = 1:numel(BIDS.subjects)
        subj = BIDS.subjects(idx).name;
        ses = BIDS.subjects(idx).session;
        mc_name = fullfile(output_dir, subj, ses, ...
            'pet', [subj '_' ses '_desc-mc_pet.nii.gz']);

        output_mean = fullfile(output_dir, subj, ses, ...
            'pet', [subj '_' ses '_desc-mc_mean.nii.gz']);

        unix(['mri_concat '  mc_name ...
            ' --mean' ...
            ' --o ' output_mean]);

        mc_mean = struct;
        mc_mean.Description = 'Mean PET image averaged across all time frames';
        mc_mean.Sources = erase(mc_name,pwd);
        mc_mean.SoftwareName = 'FreeSurfer-mri_concat';
        [~, mc_mean.SoftwareVersion] = unix('mri_concat --version');
        mc_mean.CommandLine = ['mri_concat '  erase(mc_name,pwd) ...
            ' --mean' ...
            ' --o ' erase(output_mean,pwd)];

        bids.util.jsonwrite(fullfile(output_dir, subj, ses, ...
            'pet', [subj '_' ses '_desc-mc_mean.json']),mc_mean);

        output_lta = fullfile(output_dir, subj, ses, ...
            'pet', [subj '_' ses '_from-pet_to-T1w_reg.lta']);

        unix(['mri_coreg --s '  subj ...
            ' --mov ' output_mean ...
            ' --reg ' output_lta]);

        % Write metadata to json
        coreg = struct;
        coreg.Description = 'Registration file to go from PET space to ANAT space';
        coreg.Sources = {erase(output_mean,pwd), erase(brainmask_file_bids,pwd)};
        coreg.ReferenceImage = fullfile(erase(output_dir,pwd), subj, ses, ...
            'anat', [subj '_' ses '_space-T1w_desc-brain_mask.nii.gz']);
        coreg.dof = BIDS.config.preproc.coreg.dof;
        coreg.QC = {};
        coreg.SoftwareName = 'FreeSurfer-mri_coreg';
        [~, coreg.SoftwareVersion] = unix('mri_coreg --version');
        coreg.CommandLine = ['mri_coreg --s '  subj ...
            ' --mov ' erase(output_mean,pwd) ...
            ' --reg ' erase(output_lta, pwd)];

        bids.util.jsonwrite(fullfile(output_dir, subj, ses, ...
            'pet', [subj '_' ses '_from-pet_to-T1w_reg.json']),coreg);
    end
end