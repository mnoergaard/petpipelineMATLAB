function [TimeTAC,Weights,RefTAC,RoiTAC,Name,Par,DocTxt,ImageFN]=KinmodEvaluateInpArg(KinModFunc,varargin)
%
% [TimeTAC,Weights,RefTAC,RoiTAC,Name,Time/k2ref,DocTxt]=KinmodEvaluateInpArg(KinModFunctac,varargin)
%
%  General function that can be used for evaluating the input to all the 
%  kinetic modeling functions like CreateSRTM, CreateESRTM,.... 
%  
% KinModFunc - Name of the kinetic modeling function
% varargin   - List of input to the kinetic modeling function
%
% Inputs:
% tacFN    - *.tac file from CreatePMODtac_?, contains infor about regions
%            and timing (timing should be in sec)
% Time     - structure with time information of which datapoints to use (or 
%            length of time for fitting), and eventually where challenge is
%            all of them specified in seconds
%    Time.Start - Start time for fitting
%    Time.Stop  - Stop time for fitting
%    Time.Delta - How long time used for fitting (ignored if Start 
%                 and Stop given)
%    Time.Chlng - Challenge time (only for ESRTM/ESRTM2 else ignored)
%    Time.k2ref - k2ref used for fitting (only for LoganRef else ignored).
%                 Can take values 'MRTM', 'SRTM' or a real value (iof
%                 MRTM/SRTM then these models is used for providing k2ref
%                 value [1/min].
%    Time.tstar - tstar, when to start fitting (how many samples to include) 
%                 (only for LoganRef else ignored) [min].
%
% SIFfn    - SIF file with prompts and randoms information (optional)
% Isotope  - 'C11', 'F18', 'O15', 'I123' (alternative half life of isotope
%            can be given - in minutes) (optional)
% ImageFN  - 4D image file with TAC for each voxel
%
% Par      - Normally empty, but for:
%              ESRTM - Par.Chlng
%              LoganRef - Par.k2ref, Par.tstar
%
% DocTxt   - Documentation txt to print into output files, containing info
%            about data, files used, username, server and more
%
%
% CS, 20150116
%

if ~isempty(strfind(KinModFunc,'CreateE'))
    ESRTM=1;
else
    ESRTM=0;
end
if ~isempty(strfind(KinModFunc,'LoganRef'))
    LoganRef=1;
else
    LoganRef=0;
end
if ~isempty(strfind(KinModFunc,'CreateSRTMimg'))
    EstmImage=1;
else
    EstmImage=0;
end

%
% All data has to be read from GUI
%
if nargin==1
    %
    % Read name of TAC file
    %
    [fn,pn]=uigetfile('*.tsv','Select BIDS tac file?');
    if (fn==0)
        error('%s: No PMOD.tac file selected',KinModFunc);
    end
    tacFN=[pn,fn];
    %
    % Import data
    %
    [ImpData,MinMaxTime]=ImportTACdata(tacFN);
    %
    % Decide what start/stop time for fit should be used
    %
    if ESRTM==0
        prompt={'Enter start time:','Enter stop time:'};
        name='Info about times to use';
        numlines=1;
        defaultanswer={'0',...
            num2str(MinMaxTime(2)-MinMaxTime(1))};
    else
        prompt={'Enter start time:','Enter stop time:','Enter challenge time:'};
        name='Info about times to use';
        numlines=1;
        defaultanswer={'0',...
            num2str(MinMaxTime(2)-MinMaxTime(1)),...
            num2str((MinMaxTime(2)-MinMaxTime(1))/2)};
    end
    %
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    %
    Time.Start=str2num(answer{1});
    Time.Stop=str2num(answer{2});
    if ESRTM==1
        Time.Chlng=str2num(answer{3});
    end
    %
    % Get name of SIF file, can be empty (then no weighting)
    %
    [fn,pn]=uigetfile('*.sif','Select <Noise info>.sif file?');
    if (fn==0)
        SIFfn='';
        warning('%s: No NOISE.sif file selected',KinModFunc);
    else
        SIFfn=[pn,fn];
    end
    %
    % If weighting get information on isotope
    %
    if ~isempty(SIFfn)
        sel_txt={'C11', 'F18', 'O15', 'I123'};
        [s,v] = listdlg('PromptString','What isotope?', ...
            'SelectionMode','single',...
            'ListString',sel_txt);
        if isempty(s)
            error('%s: No weighting possible without a known isotope',KinModFunc);
        else
            Isotope=sel_txt{s};
        end
    end
    %
    if LoganRef==1
        sel_txt={'MRTM', 'SRTM', 'Manual enter'};
        [s,v] = listdlg('PromptString','What method for getting k2''', ...
            'SelectionMode','single',...
            'ListString',sel_txt);
        if isempty(s)
            warning('%s: No value for k2'' entered, assumed to be huge',KinModFunc);
            Time.k2ref=1e6;
        elseif s==1||s==2
            Time.k2ref=sel_txt{s};
        else
            prompt={'Enter k2'' [1/min]:'};
            name='k2'' to use';
            numlines=1;
            defaultanswer={'0.05'};
            %
            answer=inputdlg(prompt,name,numlines,defaultanswer);
            %
            Time.k2ref=str2double(answer{1});
        end
        prompt={'Enter tstar'' [min]:'};
        name='tstar to use';
        numlines=1;
        defaultanswer={'10'};
        %
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        %
        Time.tstar=str2double(answer{1});
    end
    %
    if EstmImage==1
        [fn,pn]=uigetfile('*.img','Select <4D image>.img file?');
        if (fn==0)
            ImageFN='';
            error('%s: No ImageFN.img file selected',KinModFunc);
        else
            ImageFN=[pn,fn];
        end
    end
        
else
    %
    % Import data
    %
    tacFN=varargin{1};
    [ImpData,MinMaxTime]=ImportTACdata(tacFN);
end
%
% Evaluate or create Time argument if called with arguments
%
if nargin>1
    if nargin>=3
        if ~isempty(varargin{2})
            Time=varargin{2};
            if isfield(Time,'Start')
                if ~isfield(Time,'Stop')
                    error('%s: time argument has to have both start and stop time defined',KinModFunc);
                end
            elseif isfield(Time,'Delta')
                Time.Start=0;
                Time.Stop=0+Time.Delta;
                %Time.Chlng=(Time.Stop-Time.Start)/2;
            else
                error('%s: time argument different from blank but has no required fields defined',KinModFunc);
            end
            if ~isfield(Time,'Chlng')
                Time.Chlng=(Time.Stop-Time.Start)/2;
            end
        else
            %
            % All samples will be used
            %
            Time.Start=0;
            Time.Stop=max(ImpData.data(:,2)-ImpData.data(:,1));
            Time.Chlng=(Time.Stop-Time.Start)/2;
        end
    else
        error('KinmodEvaluateInpArg: Has to be called with either 1 or >=3 parameters - last can be empty')
    end
end
if Time.Start<0
    error('%s: Start time has to be > 0',KinModFunc);
end
Par='';
if ESRTM==1
    Par.Chlng=Time.Chlng/60;
else
    Par.Chlng=NaN;
end
if LoganRef==1
    if ~isfield(Time,'k2ref')
        Time.k2ref=1e6;
    end
    Par.k2ref=Time.k2ref;
    %
    if ~isfield(Time,'tstar')
        Time.tstar=10;
    end
    Par.tstar=Time.tstar;
else
    Par.k2ref=NaN;
    Par.tstar=nan;
end
%
%
%
if (nargin>=4)
    SIFfn=varargin{3};
end
if (nargin>=5)
    Isotope=varargin{4};
end 
if (nargin>=6)
    ImageFN=varargin{5};
end
%
% Do calculation of weights if SIFfn exist and is different from '' else
% equal weighting
%
if exist('SIFfn','var')==1 && ~isempty(SIFfn)
    %
    SIFdata=ImportSIF(SIFfn);
    %
    % Calculate weights
    % exp_corr_factor(t) = exp(ln(1/2)/Thalf * t)
    % Weights=  (L^2/T) / exp_corr_factor(t)^2
    % T is trues
    % L is frame length
    %
    Trues=SIFdata.Prompts-SIFdata.Randoms;
    %
    if ischar(Isotope)
        HalfLife=GetHalfLife(Isotope)*60;
    elseif isnumeric(Isotope)
        HalfLife=Isotope*60;
    else
        error('%s: Isotope has to be either string or numeric',KinModFunc);        
    end
    exp_corr_factor=exp(log(1/2)/HalfLife*((SIFdata.Start+SIFdata.Stop)/2));
    Weights=(SIFdata.Stop-SIFdata.Start).^2./Trues.*exp_corr_factor.^2;
    Weights=Weights/sum(Weights);
else
    %
    % Equal weighting
    %
    warning('%s: equal data weighting used',KinModFunc);
    %Weights=(ImpData.data(:,2)-ImpData.data(:,1));  % CS first proposal -
                                                     % weight with length of frame
    Weights=ones(size(ImpData.data(:,2)));
    Weights=Weights/sum(Weights);
end
%
% Do calculate mid-frame time and correct for the varying start time reported 
% by scanner due to this start on 10 kcounts/sec
%
TimeTAC=(ImpData.data(:,1)+ImpData.data(:,2))/2-ImpData.data(1,1);
%
% plot weights
%
figure
plot(TimeTAC,Weights,'*');
xlabel('Time (sec)')
ylabel('Fitting weights')
%
% Decide which of the samples should be used
%
StartInd=find(Time.Start<TimeTAC,1,'first');
if isempty(StartInd)
    StartInd=1;
end
StopInd=find(Time.Stop>TimeTAC,1,'last');
if isempty(StopInd)
    StopInd=length(TimeTAC);
end
%
% Due to the varying start time we have in the dataseries (delayed startup)
%
%
%  Do calculation of SRTM for all TACs
%
TimeTAC=TimeTAC(StartInd:StopInd)/60;   % All calculations performed in minutes
Weights=Weights(StartInd:StopInd);
%
RefTAC=ImpData.data(StartInd:StopInd,end);
RoiTAC=ImpData.data(StartInd:StopInd,3:end-1);
%
Name.Ref=ImpData.colheaders{end};
for i=3:length(ImpData.colheaders)-1
    Name.Roi{i-2}=ImpData.colheaders{i};
    RoiName{i-2}=Name.Roi{i-2};
end
[pn,fn,ext]=fileparts(tacFN);
Name.Label=fn;
%
% Documentation txt
%
DocTxt{1}=sprintf('#Date: %s\n',datestr(clock));
if ispc
    DocTxt{end+1}=sprintf('#User: %s\n',getenv('USERNAME'));
    DocTxt{end+1}=sprintf('#Host: %s\n',getenv('COMPUTERNAME'));
else
    DocTxt{end+1}=sprintf('#User: %s\n',getenv('USER'));
    DocTxt{end+1}=sprintf('#Host: %s\n',getenv('HOST'));
end
DocTxt{end+1}=sprintf('#Comp: %s\n',computer);
DocTxt{end+1}=sprintf('\n');
%
% Document which tac files used
%
DocTxt{end+1}=sprintf('Data modeled using %s based on the tac file:\n',KinModFunc(7:end));
[pn,fn,ext]=fileparts(tacFN);
if isempty(pn)
    pn=pwd;
end
tacFN=fullfile(pn,[fn ext]);
DocTxt{end+1}=sprintf('   %s\n',tacFN);
%
% Document what time interval used for fitting
% 
DocTxt{end+1}=sprintf('\nData available in TAC file [%i, %i] [sec]\n',MinMaxTime(1),MinMaxTime(2));
DocTxt{end+1}=sprintf('    Normalized times in TAC file [%i, %i] [sec]\n',0,MinMaxTime(2)-MinMaxTime(1));
DocTxt{end+1}=sprintf('    Normalized mean frame times used for fitting [%6.2f, %6.2f] [sec]\n',TimeTAC(1)*60,TimeTAC(end)*60);
if ESRTM==1
    DocTxt{end+1}=sprintf('    ESRTM: challenge at time %6.2f [sec]\n',Par.Chlng*60);
end
if LoganRef==1
    if isstr(Par.k2ref)
        DocTxt{end+1}=sprintf('    LoganRef: k2ref estimated using %s\n',Par.k2ref);
    else
        DocTxt{end+1}=sprintf('    LoganRef: k2ref set to [%6.2f] [1/min]\n',Par.k2ref);
    end
    DocTxt{end+1}=sprintf('    LoganRef: tstar set to %6.2f [min]\n',Par.tstar);
end
%DocTxt{end+1}=newline;
DocTxt{end+1}=sprintf('\n');;
%
% Document if data weighting used
%
if exist('Trues','var')
    %
    % Document if data weighting has been used
    %
    DocTxt{end+1}=sprintf('\nData has been weighted using prompts and randoms from the sif file:\n');
    [pn,fn,ext]=fileparts(SIFfn);
    if isempty(pn)
        pn=pwd;
    end
    SIFfn=fullfile(pn,[fn ext]);
    DocTxt{end+1}=sprintf('   %s\n',SIFfn);
    if ischar(Isotope)
        DocTxt{end+1}=sprintf('   Isotope: %s\n',Isotope);
    elseif isnumeric(Isotope)
        DocTxt{end+1}=sprintf('   Halflife: %s [min]\n',Isotope);
    end
end
DocTxt{end+1}=sprintf('\n\n');



    function [ImpData,MinMaxTime]=ImportTACdata(TACfn)
        %
        % Import TAC data
        %
        ImpData=importdata(TACfn);
        if isempty(ImpData)
            error('KinmodEvaluateInpArg: TAC data for modeling is empty');
        else
            MinMaxTime=[min(ImpData.data(:,1)) max(ImpData.data(:,2))];
            ImpData.data(:,1:2)=ImpData.data(:,1:2)-MinMaxTime(1);
        end
    end

end

