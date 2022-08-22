function [data]=ImportSIF(SIFfn)
%
% [data]=ImportSIF(SIFfn)
%
% Function that can import content of SIF file
%
% SIFfn - File name for SIF file (eventually it can contain path too
%
% data.Start   - Start of frame time (seconds)
% data.Stop    - Stop of frame time (seconds)
% data.Prompts - Prompts in the frame (counts)
% data.Randoms - Randoms in the frame (counts)
%
% Prompts = Trues + Randoms + Scatter but as Scatter is unknown Trues is
% normally calculated as Prompts - Randoms
%
% CS, 20140112
%
[pid]=fopen(SIFfn,'r');
if pid~=-1
    FirstLine=fgetl(pid);
    Space=strfind(FirstLine,' ');
    n=datevec(FirstLine(1:Space(2)),'dd-mm-yyyy HH:MM:SS');
    Rows=str2num(FirstLine(Space(2):Space(3)));
    Cols=str2num(FirstLine(Space(3):Space(4)));
    Vers=str2num(FirstLine(Space(4):end));
    %
    sif = fscanf(pid,'%f',[Cols Inf])';
    if size(sif,1)~=Rows
        warning('ImportSIF: No of data rows does not correspond to header');
    else
        data.Start=sif(:,1);
        data.Stop=sif(:,2);
        data.Prompts=sif(:,3);
        data.Randoms=sif(:,4);
    end
else
    warning('InportSIF: Not possible to open SIF file');
end