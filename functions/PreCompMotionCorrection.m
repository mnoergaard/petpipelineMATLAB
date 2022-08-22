function BIDS = MotionCorrection(BIDS)
%
%
%

output_dir = fullfile(BIDS.pth,'derivatives',BIDS.config.env.derivatives_dir);

for idx = 1:numel(BIDS.subjects)
    subj = BIDS.subjects(idx).name;
    ses = BIDS.subjects(idx).session;
    file = fullfile(BIDS.subjects(idx).path, 'pet', BIDS.subjects(idx).pet(1).filename);

    output_name = fullfile(output_dir, subj, ses, ...
        'pet', [subj '_' ses '_desc-mc_pet.nii.gz']);

    output_confounds = fullfile(output_dir, subj, ses, ...
        'pet', [subj '_' ses '_desc-confounds_timeseries.tsv']);

    mc_dir = dir([BIDS.pth '/derivatives/hmc_workflow/*' erase(ses,'ses-') '*' strrep(BIDS.subjects(idx).name(end-2:end),'-','_') ]);
    mc_file = fullfile(mc_dir.folder, mc_dir.name, 'concat_frames/mc.nii.gz');
    movement = fullfile(mc_dir.folder, mc_dir.name, 'plot_motion/movement.png');
    rotation = fullfile(mc_dir.folder, mc_dir.name, 'plot_motion/rotation.png');
    translation = fullfile(mc_dir.folder, mc_dir.name, 'plot_motion/translation.png');
    confounds = fullfile(mc_dir.folder, mc_dir.name, 'hmc_movement_output/hmc_confounds.tsv');
    
    mkdir(fullfile(output_dir, 'logs', subj))
    copyfile(mc_file, output_name);
    copyfile(confounds, output_confounds);
    copyfile(movement, fullfile(output_dir, 'logs', subj, [subj '_' ses '_movement.png']))
    copyfile(rotation, fullfile(output_dir, 'logs', subj, [subj '_' ses '_rotation.png']))
    copyfile(translation, fullfile(output_dir, 'logs', subj, [subj '_' ses '_translation.png']))

    mc = struct;
    mc.Description = 'Motion-corrected PET file';
    mc.Sources = erase(file,pwd);
    mc.ReferenceImage = 'Robust template using mri_robust_register';
    mc.CostFunction = 'ROB';
    mc.QC = '';
    mc.SoftwareName = 'PETPrep HMC workflow';
    mc.SoftwareVersion = 'v. 0.0.1';
    bids.util.jsonwrite(fullfile(output_dir, subj, ses, ...
        'pet', [subj '_' ses '_desc-mc_pet.json']),mc);
end

