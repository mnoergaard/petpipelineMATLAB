function [Err,Par,RoiName]=CreateSRTMimg(varargin)
%
% [Err,Par,RoiName]=CreateSRTMimg([tacFN,[StartStopTime,[SIFfn,Isotope]],ImageFN])
%
%  Function that creates a txt file that can be imported into database with
%  the SRTM
%  estimate. If missing frames in the middle of time series then mean
%  estimate for previous/following frame is used, if missing info
%  before/after first/last frame info in first/last frame is used
%
%  As High Bind TAC (used for SRTM k2' estimation is used first time series
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
% SIFfn    - SIF file with prompts and randoms information
% Isotope  - 'C11', 'F18', 'O15', 'I123' (alternative half life of isotope
%            can be given - in minutes)
% ImageFN  - Name of 4D image file
%
%
% CS, 20190208
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
[TimeTAC,Weights,RefTAC,RoiTAC,Name,dummy1,DocTxt,ImageFN]=KinmodEvaluateInpArg('CreateSRTMimg',varargin{:});
%
[img,hdr]=ReadAnalyzeImg(ImageFN);
if hdr.scale~=0
    img=img*hdr.scale;
    img=reshape(img,prod(hdr.dim(1:3)),hdr.dim(4));
end
%
% Decide energy in time series
%
TimeTACDelta=zeros(size(TimeTAC));
for i=1:length(TimeTAC)
    if i==1
        TimeTACDelta(1)=TimeTAC(1)*2;
        Tstop=TimeTAC(1)*2;
    else
        TimeTACDelta(i)=(TimeTAC(i)-Tstop)*2;
        Tstop=Tstop+TimeTACDelta(i);
    end
end
SumSqrImg=sum((img.*img.*repmat(TimeTACDelta',size(img,1),1)),2);
SumSqrRefVal=max(RefTAC.*RefTAC.*TimeTACDelta);
%
MedianImg=median(img,2);
MedianRefVal=median(RefTAC);
%
mask=(SumSqrImg>SumSqrRefVal*0.01)&(MedianImg>MedianRefVal*0.1);
MaskInd=find(mask~=0);
%
Iterations=2;
%
% Estimate parameters
%
R1=zeros(hdr.dim(1:3)');
k2=zeros(hdr.dim(1:3)');
BPnd=zeros(hdr.dim(1:3)');
CovR1=zeros(hdr.dim(1:3)');
Covk2=zeros(hdr.dim(1:3)');
CovBPnd=zeros(hdr.dim(1:3)');
ImgLogLike=zeros(hdr.dim(1:3)');
ImgMSE=zeros(hdr.dim(1:3)');
%
% Enable parallel
%
%parpool;
%
step=100;
loops=ceil(length(MaskInd)/step);
%
TC=TimeTAC;
W=Weights;
RT=RefTAC;
Iter=Iterations;
for cc=1:step
     NameTmp1.Roi{cc}=sprintf('ROI%i',cc);
end
%
Par=zeros(loops,step,3);
Cov=zeros(loops,step,3);   %Coefficient of variation
LogLikeR=zeros(loops,step);
MSER=zeros(loops,step);
%

NoInLastFrame=length(MaskInd)-(loops-1)*step;
%
parfor bb=1:loops
%for bb=1:loops
    fprintf('Estimating voxels (%i of %i)\n',bb,loops);
    if bb<loops
        imgtmp=img(MaskInd((bb-1)*step+1:bb*step),:);
    else
        imgtmp=zeros(step,size(img,2));
        imgtmp(1:NoInLastFrame,:)=img(MaskInd((bb-1)*step+1:end),:);
    end
    echo off
    [Err1,Par1,yest1,Par_cov1,LogLike1,AIC1,MSE1,FPE1]=srtm_estm_img(TC,W,RT,imgtmp',NameTmp1,Iter);
    echo on
    Par(bb,:,:)=Par1;
    for kk=1:step
        Cov(bb,kk,:)=sqrt(diag(squeeze(Par_cov1(kk,:,:))))./Par1(kk,:)'*100;
    end
    LogLikeR(bb,:)=LogLike1;
    MSER(bb,:)=MSE1;
end
%
for bb=1:loops
    if bb==loops
        NoInFrame=NoInLastFrame;
    else
        NoInFrame=step;
    end
    for ii=1:NoInFrame
        R1(MaskInd((bb-1)*step+ii))=Par(bb,ii,1);
        k2(MaskInd((bb-1)*step+ii))=Par(bb,ii,2);
        BPnd(MaskInd((bb-1)*step+ii))=Par(bb,ii,3);
        CovR1(MaskInd((bb-1)*step+ii))=Cov(bb,ii,1);
        Covk2(MaskInd((bb-1)*step+ii))=Cov(bb,ii,2);
        CovBPnd(MaskInd((bb-1)*step+ii))=Cov(bb,ii,3);
        ImgLogLike(MaskInd((bb-1)*step+ii))=LogLikeR(bb,ii);
        ImgMSE(MaskInd((bb-1)*step+ii))=MSER(bb,ii);
    end
end
%
%ErrFitInd=find(R1<0|R1>5|k2<0|k2>1|BPnd<0|BPnd>20);
%R1(ErrFitInd)=0;
%k2(ErrFitInd)=0;
%BPnd(ErrFitInd)=0;
%
%
for bb=1:3
    switch bb
        case 1
            imgo=R1;
            NameExt='_srtm_R1';
        case 2
            imgo=k2;
            NameExt='_srtm_k2';
        case 3
            imgo=BPnd;
            NameExt='_srtm_BPnd';
    end
    %
    hdro=hdr;
    hdro.name=[hdr.name NameExt];
    hdro.scale=max(imgo(:))/32767;
    hdro.dim=hdro.dim(1:3);
    %imgo=imgo/hdro.scale;
    hdro.lim=[1 0];
    hdro.pre=32;
    hdro.scale=1;
    WriteAnalyzeImg(hdro,imgo);
end

end
