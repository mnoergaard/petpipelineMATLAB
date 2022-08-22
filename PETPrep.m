function [] = PETPrep(data_dir, config)
% The PETPrep pipeline is a state-of-the-art PET preprocessing pipeline (wrapper),
% performing motion correction, co-registration, segmentation, partial
% volume correction and kinetic modeling on BIDS structured datasets
% containing at least one dynamic PET scan and one anatomical MRI scan. 
% 
% Parameters
% ----------
% data_dir : str
%   path to BIDS data directory
% config : str 
%   path to PETPrep configuration file
%
% Returns
% -------
% None
%
% Author: Martin Norgaard, Stanford University, 2022


cd(data_dir)
addpath(genpath(fullfile(data_dir,'code')));

% load BIDS data into matlab
BIDS = bids.layout(data_dir);

% load config file
BIDS.config = bids.util.jsondecode(['code/' config]);

% create derivatives directories
create_dirs_derivative(BIDS, ['code/' config])

config_num = regexp(config,'\d*','Match');
if ~isempty(config_num)
    BIDS.config.env.derivatives_dir = [BIDS.config.env.derivatives_dir config_num{1}];
end

% Recon-all FreeSurfer
ReconAll(BIDS)

% GTMSeg
GTMSeg(BIDS)

% FreeSurfer output to derivatives
ConvertFS2BIDS(BIDS)

% Motion Correction
if strcmp(BIDS.config.preproc.mc.precomp,'hmc_workflow')
    PreCompMotionCorrection(BIDS)
elseif strcmp(BIDS.config.preproc.mc.precomp,'no_hmc')
    NoMotionCorrection(BIDS)
else
    MotionCorrection(BIDS)
    PlotMotion(BIDS)
end

% Co-registration
CoReg(BIDS)
PlotCoReg(BIDS)

% GTMPVC
GTMPVC(BIDS)

% Convert PETsurfer output to BIDS
PETsurfer2TAC(BIDS)

% Kinetic Modeling
if ~strcmp(BIDS.config.preproc.pvc.pvc,'nopvc')
BIDS.config.preproc.kinmod.model = 'srtm';
KinMod(BIDS)

BIDS.config.preproc.kinmod.model = 'mrtm2';
KinMod(BIDS)

BIDS.config.preproc.kinmod.model = 'srtm2';
KinMod(BIDS)

BIDS.config.preproc.kinmod.model = 'loganref';
BIDS.config.preproc.kinmod.tstar = 2600;
KinMod(BIDS)
end

BIDS.config.preproc.pvc.pvc = 'nopvc';
BIDS.config.preproc.kinmod.model = 'srtm';
KinMod(BIDS)

BIDS.config.preproc.kinmod.model = 'mrtm2';
KinMod(BIDS)

BIDS.config.preproc.kinmod.model = 'srtm2';
KinMod(BIDS)

BIDS.config.preproc.kinmod.model = 'loganref';
BIDS.config.preproc.kinmod.tstar = 2600;
KinMod(BIDS)

% Surface-based analysis
PETVol2Surf(BIDS)

% Volume-based analysis
PETVol2Vol(BIDS)

    