function BIDS = ReconAll(BIDS)
% Interface to run FreeSurfer's recon-all for the BIDS dataset specified in
% the BIDS structure (input) obtained using bids-matlab.
%
% Parameters
% ----------
% BIDS : BIDS structure from bids-matlab including configuration file
%
% Date: 2022-08-15
% Author: Martin Norgaard

% FreeSurfer recon-all
session = bids.query(BIDS, 'sessions');
if numel(session)>=1
    session_fs = session{1};
end

participants = bids.query(BIDS, 'subjects');

output_dir = fullfile(BIDS.pth,'derivatives','freesurfer');

if BIDS.config.env.nproc > 1
    parpool('local',BIDS.config.env.nproc)
    parfor idx = 1:numel(participants)
        subj = BIDS.subjects(idx).name;
        input_data = bids.query(BIDS, 'data', 'sub', participants{idx}, 'suffix', 'T1w', 'ses',session_fs);
        mkdir(output_dir)
        setenv('SUBJECTS_DIR',output_dir)
        T2_data = bids.query(BIDS, 'data', 'sub', participants{idx}, 'suffix', 'T2w', 'ses',session_fs);
        if ~isempty(T2_data)
            cmd = ['recon-all -subjid sub-' num2str(participants{idx}) ...
                ' -i ' input_data{1} ...
                ' -T2 ' T2_data{1} ...
                ' -T2pial' ...
                ' -all -cw256'];
        else
            cmd = ['recon-all -subjid sub-' num2str(participants{idx}) ...
                ' -i ' input_data{1} ...
                ' -all -cw256'];
        end

        recon = struct;
        recon.Description = '';
        recon.Sources = {erase(input_data{1},pwd)};
        if ~isempty(T2_data)
            recon.Sources{end+1} = erase(T2_data{1},pwd);
        end
        recon.SoftwareName = 'FreeSurfer-recon-all';
        [~, recon.SoftwareVersion] = unix('recon-all --version');
        if ~isempty(T2_data)
            recon.CommandLine = ['recon-all -subjid sub-' num2str(participants{idx}) ...
                ' -i ' erase(input_data{1},pwd) ...
                ' -T2 ' erase(T2_data{1},pwd) ...
                ' -T2pial' ...
                ' -all -cw256'];
        else
        recon.CommandLine = ['recon-all -subjid sub-' num2str(participants{idx}) ...
            ' -i ' erase(input_data{1},pwd) ...
            ' -all -cw256'];
        end

        bids.util.jsonwrite(fullfile(BIDS.pth,'derivatives',BIDS.config.env.derivatives_dir, subj, ['ses-' session_fs], ...
            'anat', [subj '_' session_fs '_space-T1w_anat.json']),recon);

        unix(cmd)
    end
    delete(gcp('nocreate'));
else
    for idx = 1:numel(participants)
        subj = BIDS.subjects(idx).name;
        input_data = bids.query(BIDS, 'data', 'sub', participants{idx}, 'suffix', 'T1w', 'ses',session_fs);
        mkdir(output_dir)
        setenv('SUBJECTS_DIR',output_dir)
        T2_data = bids.query(BIDS, 'data', 'sub', participants{idx}, 'suffix', 'T2w', 'ses',session_fs);

        if ~isempty(T2_data)
            cmd = ['recon-all -subjid sub-' num2str(participants{idx}) ...
                ' -i ' input_data{1} ...
                ' -T2 ' T2_data{1} ...
                ' -T2pial' ...
                ' -all -cw256'];
        else
            cmd = ['recon-all -subjid sub-' num2str(participants{idx}) ...
                ' -i ' input_data{1} ...
                ' -all -cw256'];
        end

        recon = struct;
        recon.Description = '';
        recon.Sources = {erase(input_data{1},pwd)};
        if ~isempty(T2_data)
            recon.Sources{end+1} = erase(T2_data{1},pwd);
        end
        recon.SoftwareName = 'FreeSurfer-recon-all';
        [~, recon.SoftwareVersion] = unix('recon-all --version');
        if ~isempty(T2_data)
            recon.CommandLine = ['recon-all -subjid sub-' num2str(participants{idx}) ...
                ' -i ' erase(input_data{1},pwd) ...
                ' -T2 ' erase(T2_data{1},pwd) ...
                ' -T2pial' ...
                ' -all -cw256'];
        else
        recon.CommandLine = ['recon-all -subjid sub-' num2str(participants{idx}) ...
            ' -i ' erase(input_data{1},pwd) ...
            ' -all -cw256'];
        end

        bids.util.jsonwrite(fullfile(BIDS.pth,'derivatives',BIDS.config.env.derivatives_dir, subj, ['ses-' session_fs], ...
            'anat', [subj '_' session_fs '_space-T1w_anat.json']),recon);

        unix(cmd)
    end
end