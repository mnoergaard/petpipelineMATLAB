function [Err,Par,RoiName]=CreateLoganRef(varargin)
%
% [Err,Par,RoiName]=CreateLoganRef([tacFN,[StartStopTime,[SIFfn,Isotope]]])
%
%  Function that creates a txt file that can be imported into database with
%  the Logan Ref model estimate, The logan method is implemented with a
%  reference region, like in the Logan paper 1996, JCBFM 16:834-40
%   
%  If missing frames in the middle of time series then mean
%  estimate for previous/following frame is used, if missing info
%  before/after first/last frame info in first/last frame is used
%
%  As High Bind TAC (used for Logan k2' estimation is used first time series
%  (column 3 in tac file)) and as Ref TAC is used last column in tac file
% 
%  vargin is evaluated as:
% 
% tacFN    - *.tac file from CreatePMODtac_?, contains infor about regions
%            and timing (timing should be in sec)
% StartStopChlgTime     - structure with time information of which datapoints to use (or 
%            length of time for fitting), and eventually where challenge is
%            all of them specified in seconds. If empty then they are
%            attumaticllay decided from tac file, challenge in the middle
%    Time.Start - Start time for fitting
%    Time.Stop  - Stop time for fitting
%    Time.Delta - How long time used for fitting (ignored if Start 
%                 and Stop given)
%    Time.k2ref - The k2' used for estimating the BPnd, can be set to
%                 either the text strings 'MRTM' or 'SRTM' and then these
%                 models will be used for estimation of k2', or a value and
%                 then the value will be used. If not provided the value
%                 10000 will be sued (the same as ignoring this part of the
%                 LoganRef model.
%
% SIFfn    - SIF file with prompts and randoms information
% Isotope  - 'C11', 'F18', 'O15', 'I123' (alternative half life of isotope
%            can be given - in minutes)
%
%
% CS, 20180509
%

%
% Redo plot with selected regions
%
if (nargin==1)&&strcmp(varargin{1},'ShowROIs')
    ShowROIs;
    return
end
%
% Evaluate input parameters
%
[TimeTAC,Weights,RefTAC,RoiTAC,Name,ExtPar,DocTxt]=KinmodEvaluateInpArg('CreateLoganRef',varargin{:});
RoiName=Name.Roi;
%
% Estimate parameters
%
if isstr(ExtPar.k2ref)
    if strfind(ExtPar.k2ref,'MRTM')
        [Err,Par,yest,Par_cov,LogLike,AIC, MSE, FPE]=mrtm_estm(TimeTAC,Weights,RefTAC,RoiTAC,Name);
        ExtPar.k2ref=Par{1}(2);
        ExtPar.k2ref_cov=Par_cov{1}(2,2);
    else
        [Err,Par,yest,Par_cov,LogLike,AIC, MSE, FPE]=srtm_estm(TimeTAC,Weights,RefTAC,RoiTAC,Name);
        ExtPar.k2ref=Par{1}(2)/Par{1}(1);
        ExtPar.k2ref_cov=Par_cov{1}(2,2)/(Par{1}(1,1)^2)+Par_cov{1}(1,1)*(Par{1}(2)^2)/(Par{1}(1)^4);
    end
else
    ExtPar.k2ref_cov=nan;
end
%
fprintf('Used k2ref is: %6.4f, and tstar: %6.2f\n',ExtPar.k2ref,ExtPar.tstar);
[Err,Par,yest,Par_cov,LogLike,AIC,MSE,FPE]=logan_ref_estm(TimeTAC,Weights,RefTAC,RoiTAC,Name,ExtPar);
%
% Visualize and save output
%
%DataTAC.Time=TimeTAC;
%DataTAC.Ref=RefTAC;
%DataTAC.Roi=RoiTAC;
%DataTAC.RoiEst=yest;
%MakePlot(DataTAC,Name,-1);
%
% Print data to screen and file
%
fprintf('\n\nData for Logan ref model\n\n');
fprintf('   Data from %s\n\n',Name.Label');
for i=1:size(RoiTAC,2)
    fprintf('   Logan Ref: %s: Kappa2=%5.3f, BPnd=%5.3f, k2ref=%5.3f, WeightErr: %6.2e, MSE: %6.2e, FPE: %6.2e, LogLikelihood: %6.3e, AIC: %6.2e\n',...
        Name.Roi{i},Par{i},ExtPar.k2ref,Err(i),MSE(i),FPE(i),LogLike(i),AIC(i));
end
%
[pn,fn,ext]=fileparts(Name.Label);
print('-djpeg','-r300',[fn(1:end-4) 'km-loganref']);
%
% Save data to file
%
[pn,fn,ext]=fileparts(Name.Label);
pid=fopen([fn(1:end-4) 'km-loganref.tsv'],'w');
if pid==-1
    warning('Logan: not possible to open file for saving Logan reference results - probably write protected directory\n\n');
else
    %
    % Document date user and so on
    %
    %for i=1:length(DocTxt)
    %    fprintf(pid,'%s',DocTxt{i});
    %end
    %
    % Write data
    %
    fprintf(pid,'Logan Ref\tKappa2\tBPnd\tk2ref\tMSE\tLogLikelihood\tAIC\t%%Cov Kappa2\t%%Cov BPnd\t%%Cov k2ref\n');
    for i=1:size(RoiTAC,2)
        fprintf(pid,'%s\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n',...
            Name.Roi{i},Par{i},ExtPar.k2ref,MSE(i),LogLike(i),AIC(i),sqrt(diag(Par_cov{i}))./Par{i}*100,ExtPar.k2ref_cov/ExtPar.k2ref*100);
    end
    fclose(pid);
end


    function MakePlot(DataTAC,Name,ROIs)
        %
        fig_h=findobj('type','figure','tag','CreateLoganfig');
        if ROIs==-1
            New=1;
            delete(fig_h);
        else
            New=0;
        end
        if (nargin==2)||New==1
            ROIs=1:size(DataTAC.Roi,2);
        end
        %
        if isempty(fig_h)||New==1
            h=figure;
            set(h,'tag','CreateLoganfig');
            ud.Data=DataTAC;
            ud.Name=Name;
            set(h,'UserData',ud);
            mh=uimenu(h,'Label','Logan');
            seh1=uimenu(mh,'Label','Select ROIs','Accelerator','S','CallBack','CreateLogan(''ShowROIs'');');
        else
            ud=get(fig_h,'UserData');
            DataTAC=ud.Data;
            Name=ud.Name;
        end
        subplot(211),plot(DataTAC.Time,DataTAC.Roi(:,ROIs),'*',DataTAC.Time,DataTAC.Ref,'+');
        hold on
        d=ver('matlab');
        if (d.Version(1)>'8')
            set(gca,'ColorOrderIndex',1);
        end
        plot(DataTAC.Time,[DataTAC.RoiEst(:,ROIs), DataTAC.Ref]);
        grid
        hold off
        xlabel('Time [min]');
        ylabel('Activity [Bq/cc]')
        h=title(sprintf('Logan fit: %s (%s)' ,Name.Label,datestr(clock,0)));
        set(h,'interpreter','none');
        ha1=gca;
        %
        subplot(212),plot(DataTAC.Time,DataTAC.Roi(:,ROIs)-DataTAC.RoiEst(:,ROIs));
        xlabel('Time [min]');
        ylabel('Residuals: Activity [Bq/cc]')
        ha2=gca;
        grid
        %
        if length(ROIs)<16    % More than 16 ROI's then not possible to do legend
            nn=Name.Roi(ROIs);
            nn{end+1}=Name.Ref;
            hl=legend(ha1,nn,'Location','EastOutside');
            h=findobj('type','text');
            set(h,'fontsize',8);
            set(h,'interpreter','none');
            %
            %newPosition = [0.4 0.4 0.2 0.2];
            %newUnits = 'normalized';
            %set(hl,'Position', newPosition,'Units', newUnits);
            pos1=get(ha1,'Position');
            pos2=get(ha2,'Position');
            pos2(3)=pos1(3);
            %set(ha1,'Position',pos1);
            set(ha2,'Position',pos2);
            %
            pos_l=get(hl,'Position');
            pos_l(3)=pos_l(3)*1.2;
            pos_l(2)=0.5-0.5*pos_l(4);
            set(hl,'Position',pos_l);
            %
            set(ha1,'Position',pos1);
        end
        
        h=findobj('type','text');
        set(h,'interpreter','none');
    end

    function ShowROIs
        fig_h=findobj('type','figure','tag','CreateLoganfig');
        ud=get(fig_h,'UserData');
        DataTAC=ud.Data;
        Name=ud.Name;
        %
        [s,v] = listdlg('PromptString','Select ROIs:',...
            'SelectionMode','multiple',...
            'ListString',Name.Roi);
        %
        MakePlot(DataTAC,Name,s);
    end
end
