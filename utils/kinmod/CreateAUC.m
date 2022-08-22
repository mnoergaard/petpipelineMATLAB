function Problem=CreateAUC(tacFN,StartStopTime)
%
% Err=CreateAUC([tacFN,[StartStopTime]])
%
%  Function that creates a txt file that can be imported into database with the AUC 
%  estimate. If missing frames in the middle of time series then mean
%  estimate for previous/following frame is used, if missing info
%  before/after first/last frame info in first/last frame is used
% 
% tacFN    - *.tac file from CreatePMODtac_?, contains infor about regions
%            and timing
% StartStopTime - two values defining from where to when at time axis to do
%            integration. If only 1 number giver this is assumed to be the
%            length of the integration from first time point (like 2 hours for SB acquisition) 
%
% Err      - Error (0/1) if missing data start/end/middle set to 1
%
% CS, 20131129
%
Problem=0;
if nargin==0
    [fn,pn]=uigetfile('*.tac','Select <PMOD>.tac file?');
    if (fn==0)
        error('CreateAUC: No PMOD.tac file selected');
    end
    tacFN=[pn,fn];
end
%
ImpData=importdata(tacFN);
%
if nargin==2
    if ~isnumeric(StartStopTime)
        error('CreateAUC: Second argument has to be a 2 value argument');
    else
        if length(StartStopTime)==1
            StartStopTime=[min(ImpData.data(:,1)) min(ImpData.data(:,1))+StartStopTime];
        elseif length(StartStopTime)>2
            error('CreateAUC: Second argument has to be a 2 value argument');
        end
    end
else
    StartStopTime=[min(ImpData.data(:,1)) max(ImpData.data(:,2))];
end
if StartStopTime(1)<0
    error('Start time for AUC cal has to be > 0');
end
%
% Due to the varying start time we have in the dataseries (delayed startup)
% we subtract 
%
%
%  Do calculation of AUC for all TACs
%
IntTime=StartStopTime(2)-StartStopTime(1);
DeltaTime=zeros(size(ImpData.data,1),1);
for i=1:length(DeltaTime)
    if (min(StartStopTime)<=ImpData.data(i,1))
        if (max(StartStopTime)>=ImpData.data(i,2))
            DeltaTime(i)=ImpData.data(i,2)-ImpData.data(i,1);
        else
            DeltaTime(i)=max(StartStopTime)-ImpData.data(i,1);
        end
    else
        if (max(StartStopTime)>=ImpData.data(i,2))
            DeltaTime(i)=ImpData.data(i,2)-min(StartStopTime);
        else
            DeltaTime(i)=max(StartStopTime)-min(StartStopTime);
        end
    end
end
DeltaTime(DeltaTime<0)=0;
%
IntData=zeros(1,size(ImpData.data,2)-2);
for i=3:size(ImpData.data,2)
    IntData(i-2)=sum(ImpData.data(:,i).*DeltaTime);
end
%
% Looking for missing data in the middle of time series and interpolating
% these and adding them to integral
% Stop time for previous frame should be start time for next frame
%
CompTime=ImpData.data(1:end-1,2)==ImpData.data(2:end,1);
pos=find(CompTime==0);
for i=1:length(pos)
    DTime=ImpData.data(pos(i)+1,1)-ImpData.data(pos(i),2);
    Val=(ImpData.data(pos(i),3:end)+ImpData.data(pos(i)+1,3:end))/2;
    DVal=Val*DTime;
    IntData=IntData+DVal;
    Problem=1;
    fprintf('CreateAUC: missing data in middle of time series for: %s\n',tacFN);
end
%
% Looking if integration interval is before/after last frame, and
% exentually extrapolation these and integrating it up
% For start time it assumes that value at time 0 is 0 for all ROI's
%
if (StartStopTime(1)<min(ImpData.data(:,1)))
    DTime=StartStopTime(1)-min(ImpData.data(:,1));
    Val=(ImpData.data(1,3:end)+0)/2;    %Assume start in 0,0
    DVal=Val*DTime;
    IntData=IntData+DVal;
    Problem=1;
    fprintf('CreateAUC: missing data in start of time series for: %s\n',tacFN);
end
if (StartStopTime(2)>max(ImpData.data(:,2)))
    DTime=StartStopTime(2)-max(ImpData.data(:,2));
    Val=ImpData.data(end,3:end);
    DVal=Val*DTime;
    IntData=IntData+DVal;
    Problem=1;
    fprintf('CreateAUC: missing data in end of time series for: %s\n',tacFN);
end
%
IntData=IntData/IntTime;
%
[pn,fn,ext]=fileparts(tacFN);
pid=fopen(fullfile(pn,[fn '.auc']),'w');
if pid==-1
    error('Not able to create file with AUC values');
else
    fprintf(pid,'ROI\tnormAUC\tStart time\tStop time\tProblem\n');
    for i=1:length(IntData)
        fprintf(pid,'%s\t%e\t%e\t%e\t%i\n',ImpData.colheaders{i+2},IntData(i),StartStopTime,Problem);
    end
end
fclose(pid);
