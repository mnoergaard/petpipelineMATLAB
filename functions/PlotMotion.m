function BIDS = PlotMotion(BIDS)
%
%
%

output_dir = fullfile(BIDS.pth,'derivatives',BIDS.config.env.derivatives_dir);

for idx = 1:numel(BIDS.subjects)
    subj = BIDS.subjects(idx).name;
    ses = BIDS.subjects(idx).session;
    file = fullfile(BIDS.subjects(idx).path, 'pet', BIDS.subjects(idx).pet(1).filename);

    logs_dir = fullfile(output_dir, 'logs', subj);

    mkdir(logs_dir);

    mc_plot_file = fullfile(output_dir, subj, ses, ...
        'pet', [subj '_' ses '_desc-mc_pet.nii.gz.par']);

    mc_data = load(mc_plot_file);
    mc.rot_x = mc_data(:,1); 
    mc.rot_y = mc_data(:,2); 
    mc.rot_z = mc_data(:,3); 
    mc.trans_x = mc_data(:,4); 
    mc.trans_y = mc_data(:,5); 
    mc.trans_z = mc_data(:,6);

    metadata = bids.query(BIDS, 'metadata', 'sub', subj, 'ses', ses, ...
        'modality', 'pet', 'suffix', 'pet');

    time = cumsum(metadata.FrameDuration);

    figure; 
    subplot(2,1,1)
    plot(time,mc_data(:,1:3),'.-'); grid on; hold on;
    plot(time, zeros(size(time)),'k-');
    xlim([0 max(time)])
    xlabel('Time [sec]')
    ylabel('Rotation [degrees]');
    legend('X','Y','Z');
    subplot(2,1,2)
    plot(time,mc_data(:,4:6),'.-'); grid on; hold on;
    plot(time, zeros(size(time)),'k-')
    xlim([0 max(time)])
    xlabel('Time [sec]')
    ylabel('Translation [mm]');
    legend('X','Y','Z');
    
    % save image
    saveas(gcf, fullfile(logs_dir,[subj '_' ses '_motion.png']));
    close(gcf)

    % save confounds
    bids.util.tsvwrite(fullfile(output_dir, subj, ses, ...
        'pet', [subj '_' ses '_desc-confounds_timeseries.tsv']),mc);

    % update corresponding json file to include QC
    mc_meta = bids.util.jsondecode(fullfile(output_dir, subj, ses, ...
                'pet', [subj '_' ses '_desc-mc_pet.json']));

    mc_meta.QC = fullfile(erase(logs_dir,pwd),[subj '_' ses '_motion.png']);

    bids.util.jsonwrite(fullfile(output_dir, subj, ses, ...
                'pet', [subj '_' ses '_desc-mc_pet.json']),mc_meta);

end

