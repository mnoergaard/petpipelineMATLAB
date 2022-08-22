function BIDS = SmoothPET(BIDS)
%
%
%

output_dir = fullfile(BIDS.pth,'derivatives',BIDS.config.env.derivatives_dir);
setenv('SUBJECTS_DIR',fullfile(BIDS.pth,'derivatives','freesurfer'));

