function BIDS = MotionCorrection(BIDS)
%
%
%

output_dir = fullfile(BIDS.pth,'derivatives',BIDS.config.env.derivatives_dir);

for idx = 1:numel(BIDS.subjects)
    subj = BIDS.subjects(idx).name;
    ses = BIDS.subjects(idx).session;

    if ~isempty(BIDS.subjects(idx).pet)
    file = fullfile(BIDS.subjects(idx).path, 'pet', BIDS.subjects(idx).pet(1).filename);

    output_name = fullfile(output_dir, subj, ses, ...
        'pet', [subj '_' ses '_desc-mc_pet.nii.gz']);

    output_json = fullfile(output_dir, subj, ses, ...
        'pet', [subj '_' ses '_desc-mc_pet.json']);

    output_confounds = fullfile(output_dir, subj, ses, ...
        'pet', [subj '_' ses '_desc-confounds_timeseries.tsv']);


    mc_file = fullfile([BIDS.pth '/derivatives/petprep_hmc/' subj '/' ses '/' subj '_' ses '_desc-mc_pet.nii.gz']);
    mc_json = fullfile([BIDS.pth '/derivatives/petprep_hmc/' subj '/' ses '/' subj '_' ses '_desc-mc_pet.json']);
    movement = fullfile([BIDS.pth '/derivatives/petprep_hmc/' subj '/' ses '/' subj '_' ses '_movement.png']);
    rotation = fullfile([BIDS.pth '/derivatives/petprep_hmc/' subj '/' ses '/' subj '_' ses '_rotation.png']);
    translation = fullfile([BIDS.pth '/derivatives/petprep_hmc/' subj '/' ses '/' subj '_' ses '_translation.png']);
    confounds = fullfile([BIDS.pth '/derivatives/petprep_hmc/' subj '/' ses '/' subj '_' ses '_desc-confounds_timeseries.tsv']);
    nmc_plot = fullfile([BIDS.pth '/derivatives/petprep_hmc/' subj '/' ses '/' subj '_' ses '_desc-without_motion_correction.gif']);
    mc_plot = fullfile([BIDS.pth '/derivatives/petprep_hmc/' subj '/' ses '/' subj '_' ses '_desc-with_motion_correction.gif']);

    if exist(mc_file,"file")
    mkdir(fullfile(output_dir, 'logs', subj))
    copyfile(mc_file, output_name);
    copyfile(mc_json, output_json);
    copyfile(confounds, output_confounds);
    copyfile(movement, fullfile(output_dir, 'logs', subj, [subj '_' ses '_movement.png']))
    copyfile(rotation, fullfile(output_dir, 'logs', subj, [subj '_' ses '_rotation.png']))
    copyfile(translation, fullfile(output_dir, 'logs', subj, [subj '_' ses '_translation.png']))
    copyfile(nmc_plot, fullfile(output_dir, 'logs', subj, [subj '_' ses '_desc-without_motion_correction.gif']))
    copyfile(mc_plot, fullfile(output_dir, 'logs', subj, [subj '_' ses '_desc-with_motion_correction.gif']))

    end
    end
end

