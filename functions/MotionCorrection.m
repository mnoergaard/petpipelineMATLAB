function BIDS = MotionCorrection(BIDS)
%
%
%

output_dir = fullfile(BIDS.pth,'derivatives',BIDS.config.env.derivatives_dir);

if BIDS.config.env.nproc > 1
    parpool('local',BIDS.config.env.nproc)
    parfor idx = 1:numel(BIDS.subjects)
            subj = BIDS.subjects(idx).name;
            ses = BIDS.subjects(idx).session;
            file = fullfile(BIDS.subjects(idx).path, 'pet', BIDS.subjects(idx).pet(1).filename);

            output_name = fullfile(output_dir, subj, ses, ...
                'pet', [subj '_' ses '_desc-mc_pet.nii.gz']);

            unix(['mcflirt -in '  file ...
                ' -out ' output_name ...
                ' -cost ' BIDS.config.preproc.mc.cost ...
                ' -refvol ' num2str(BIDS.config.preproc.mc.refvol) ...
                ' -plots']);
            
            mc = struct;
            mc.Description = 'Motion-corrected PET file';
            mc.Sources = erase(file,pwd);
            mc.ReferenceImage = [num2str(BIDS.config.preproc.mc.refvol)];
            mc.CostFunction = BIDS.config.preproc.mc.cost;
            mc.QC = '';
            mc.SoftwareName = 'FSL-mcflirt';
            mc.SoftwareVersion = 'v. 6.0';
            mc.CommandLine = ['mcflirt -in '  erase(file,pwd) ...
                ' -out ' erase(output_name, pwd) ...
                ' -cost ' BIDS.config.preproc.mc.cost ...
                ' -refvol ' num2str(BIDS.config.preproc.mc.refvol) ...
                ' -plots'];
            bids.util.jsonwrite(fullfile(output_dir, subj, ses, ...
                'pet', [subj '_' ses '_desc-mc_pet.json']),mc);
    end
    delete(gcp('nocreate'));
else
    for idx = 1:numel(BIDS.subjects)
            subj = BIDS.subjects(idx).name;
            ses = BIDS.subjects(idx).session;
            file = fullfile(BIDS.subjects(idx).path, 'pet', BIDS.subjects(idx).pet(1).filename);

            output_name = fullfile(output_dir, subj, ses, ...
                'pet', [subj '_' ses '_desc-mc_pet.nii.gz']);

            unix(['mcflirt -in '  file ...
                ' -out ' output_name ...
                ' -cost ' BIDS.config.preproc.mc.cost ...
                ' -refvol ' num2str(BIDS.config.preproc.mc.refvol) ...
                ' -plots']);
            
            mc = struct;
            mc.Description = 'Motion-corrected PET file';
            mc.Sources = erase(file,pwd);
            mc.ReferenceImage = [num2str(BIDS.config.preproc.mc.refvol)];
            mc.CostFunction = BIDS.config.preproc.mc.cost;
            mc.QC = '';
            mc.SoftwareName = 'FSL-mcflirt';
            mc.SoftwareVersion = 'v. 6.0';
            mc.CommandLine = ['mcflirt -in '  erase(file,pwd) ...
                ' -out ' erase(output_name, pwd) ...
                ' -cost ' BIDS.config.preproc.mc.cost ...
                ' -refvol ' num2str(BIDS.config.preproc.mc.refvol) ...
                ' -plots'];
            bids.util.jsonwrite(fullfile(output_dir, subj, ses, ...
                'pet', [subj '_' ses '_desc-mc_pet.json']),mc);
    end
end

