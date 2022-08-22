function BIDS = KinMod(BIDS)
%
%
%

output_dir = fullfile(BIDS.pth,'derivatives',BIDS.config.env.derivatives_dir);

if BIDS.config.env.nproc > 1
    parpool('local',BIDS.config.env.nproc)
    parfor idx = 1:numel(BIDS.subjects)
        subj = BIDS.subjects(idx).name;
        ses = BIDS.subjects(idx).session;
        tac_file = fullfile(output_dir, subj, ses, ...
        'pet', [subj '_' ses '_pvc-' ...
        BIDS.config.preproc.pvc.pvc ...
        '_desc-mc_tacs.tsv']);

        Time = struct('Start', BIDS.config.preproc.kinmod.start_time, ...
            'Stop', BIDS.config.preproc.kinmod.end_time);

        model_output_dir = fullfile(output_dir, subj, ses, ...
            'pet');

        cd(model_output_dir)
        switch BIDS.config.preproc.kinmod.model

            case 'mrtm2'
                CreateMRTM2(tac_file, Time);

            case 'srtm'
                CreateSRTM(tac_file, Time);

            case 'srtm2'
                CreateSRTM2(tac_file, Time);

            case 'loganref'
                Time.tstar = BIDS.config.preproc.kinmod.tstar / 60;
                CreateLoganRef(tac_file, Time);

        end
        close all
        cd(BIDS.pth)

    end
    delete(gcp('nocreate'));
else
    for idx = 1:numel(BIDS.subjects)
        subj = BIDS.subjects(idx).name;
        ses = BIDS.subjects(idx).session;
        tac_file = fullfile(output_dir, subj, ses, ...
        'pet', [subj '_' ses '_pvc-' ...
        BIDS.config.preproc.pvc.pvc ...
        '_desc-mc_tacs.tsv']);

        tacs = bids.util.tsvread(tac_file);
        Time = struct('Start', tacs.frame_start, 'Stop', tacs.frame_end);

        model_output_dir = fullfile(output_dir, subj, ses, ...
            'pet');

        cd(model_output_dir)
        switch BIDS.config.preproc.kinmod.model

            case 'mrtm2'
                CreateMRTM2(tac_file, Time);

            case 'srtm'
                CreateSRTM(tac_file, Time);

            case 'srtm2'
                CreateSRTM2(tac_file, Time);

            case 'loganref'
                Time.tstar = BIDS.config.preproc.kinmod.tstar / 60;
                CreateLoganRef(tac_file, Time);

        end
        close all
        cd(BIDS.pth)

    end
end
