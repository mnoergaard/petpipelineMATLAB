function [Err,Par,yest,Par_cov,LogLike,AIC, MSE, FPE]=mrtm2_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name,k2p)
%
%  function  [Err,Par,yest]=mrtm_estm(TimeTAC,FrameWeight,RefTAC,RoiTACName,k2p)
%
% Implementation of MRTM estimation
%
%   TimeTAC   - mean time for each frame (in min, assumes that frame 1 starts at
%               time zero, so e.g. if length of first frame is 1/6 min, then
%               first time point should be 1/12 min)
%   FrameWeight - Weight of frame when fitting:
%                   W=(FrameLength^2/Trues)*exp_correct_factor
%   RefTAC    - Reference tissue TAC
%   RoiTAC    - Roi tissue TAC (can contain multiple tissue curves as columns)
%   Name      - Name.Ref - String with name of ref roi
%               Name.Roi - Cell array with names of regions
%               Name.Label - Label string for data series 
%   k2p       - k2' used when estimating BPnd and R1 for all regions (common k2 for the reference region)
%
%   Par -      [R1, BPnd]
%
% CS, 20140805
%
%
xdata=[TimeTAC RefTAC FrameWeight];
%
%
for i=1:size(RoiTAC,2)
    ydata=RoiTAC(:,i);
    if (size(RoiTAC,2) < 100) ||...    % Regional analysis
       (rem(i,10000) == 0)             % Voxel based analysis 
        fprintf('\nEstimating ROI (%i/%i): %s\n',i,size(RoiTAC,2),Name.Roi{i});
    end
    [Err(i), Par{i}, yest(:,i), Par_cov{i}, LogLike(i), AIC(i), MSE(i), FPE(i)]=EstmMRTM2Tissue;
end
%

    function [Err, Par, yest, Par_cov, LogLike, AIC, MSE, FPE]=EstmMRTM2Tissue
        %
        %
        % From kinetic course notes
        %
        % Ct(T)=R1 * K2' * [Int[0;T] C'(t)dt + 1/k2'*C'(T)] - k2 * Int[0;T] C(t)dt
        % BPnd=-(b1/b2+1); R1=b1/k2';
        %
        Crefmod_int=KinmodCumtrapz_l(xdata(:,1),RefTAC)+1/k2p*RefTAC;
        Croi_int=KinmodCumtrapz_l(xdata(:,1),ydata);
        %
        x=[Crefmod_int,Croi_int];
        y=ydata;
        %
        % Linear glm model
        %
        [b,dev,stats] = glmfit(x, y, 'normal', 'constant', 'off', 'weights',xdata(:,3));
        lin.yhat=x*b;
        lin.Err=sum((y-lin.yhat).^2.*xdata(:,3))/(length(y)-2);
        lin.R1=b(1)/k2p;
        lin.BPnd=-(b(1)/b(2)+1);
        lin.R1_var=stats.covb(1,1)/k2p^2;
        lin.BPnd_var=(stats.covb(1,1)*b(2)^2-2*stats.covb(1,2)*b(1)*b(2)+stats.covb(2,2)*b(1)^2)/...
            b(2)^4;
        Par=[lin.R1 lin.BPnd];
        Par_cov=zeros(2);
        Par_cov(1,1)=lin.R1_var;
        Par_cov(2,2)=lin.BPnd_var;
        Err=lin.Err;
        yest=lin.yhat;
        %
        % Calc stats
        %
        MSE=sum((yest-ydata).^2)/(length(ydata)-4); % 3 par + std err
        FPE=sum((yest-ydata).^2)*(length(ydata)+4)/(length(ydata)-4); % 3 par + std err
        SigmaSqr=std(yest-ydata)^2;
        LogLike=-0.5*length(ydata)*log(2*pi*SigmaSqr)-0.5*sum((yest-ydata).^2)/SigmaSqr;
        AIC=-2*LogLike+2*4;  % From Klaus Holst. 4 parameters is 3 model parameters (2 from SRTM2 + k2' from SRTM) and noise variance
    end


    function Ct=MRTM2(p,xd)
        %
        % Model of Reference tissue solution  (eq 6)
        %
        R1=p(1);
        BPnd=p(2);
        %
        % Ct(T)=R1 * K2' * [Int[0;T] C'(t)dt + 1/k2'*C'(T)] - k2 * Int[0;T] C(t)dt
        % BPnd=-(b1/b2+1); R1=b1/k2';
        %
        % b1=R1*k2'; b2=-b1/(1+BPnd)=-(R1*k2p)/(BPnd+1)
        %
        Ct=R1*k2p*xd(:,1) - (R1*k2p)/(1+BPnd)*xd(:,2);
        %

    end        
        
end



