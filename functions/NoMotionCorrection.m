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
    
    mkdir(fullfile(output_dir, 'logs', subj))
    copyfile(file, output_name);

    mc = struct;
    mc.Description = 'Non motion-corrected PET file';
    mc.Sources = erase(file,pwd);
    mc.ReferenceImage = '';
    mc.CostFunction = '';
    mc.QC = '';
    mc.SoftwareName = '';
    mc.SoftwareVersion = '';
    bids.util.jsonwrite(fullfile(output_dir, subj, ses, ...
        'pet', [subj '_' ses '_desc-mc_pet.json']),mc);
end

