function [Err,Par,yest,Par_cov,LogLike,AIC, MSE, FPE]=srtm2_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name,k2p)
%
%  function  [Err,Par,yest]=srtm_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name[,Iter])
%
% Implementation of SRTM estimation
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
% CS, 20140825
%
%
Iter=10;
%
xdata=[TimeTAC RefTAC FrameWeight];
%
for i=1:size(RoiTAC,2)
    ydata=RoiTAC(:,i);
    fprintf('\nEstimating ROI (%i/%i): %s\n',i,size(RoiTAC,2),Name.Roi{i});
    [Err(i), Par{i}, yest(:,i), Par_cov{i}, LogLike(i), AIC(i), MSE(i), FPE(i)]=EstmSRTM2Tissue;
end
%

    function [Err, Par, yest, Par_cov, LogLike, AIC, MSE, FPE]=EstmSRTM2Tissue
        %
        %
        for k=1:Iter
            fprintf('  Iteration: %i',k);
            Par0=(1+0.1*randn(1,2)).*[1.0,2];   % R1, BP0
            optims.MaxIter=1000;
            %[Par_n{k}]=fminsearch(@CostRefTissue,Par0,optims);
            %[Err_n(k),yest_n{k}]=CostRefTissue(Par_n{k});
            opt=statset;
            opt.FunValCheck='off';
            %[Par_n{k},R,J,Par_cov_n{k},Err_n(k)]=nlinfit(xdata(:,1:2),ydata,@SRTM2,Par0,opt);
            [Par_n{k},R,J,Par_cov_n{k},Err_n(k)]=nlinfit(xdata(:,1:2),ydata,@SRTM2,Par0,opt,'Weights',xdata(:,3));
            yest_n{k}=SRTM2(Par_n{k},xdata(:,1:2));
            fprintf(', LSQ: %6.3e\n',Err_n(k));
        end
        [minv,ind]=min(Err_n);
        Par=[Par_n{ind}];
        Par_cov=[Par_cov_n{ind}];
        Err=Err_n(ind);
        yest=yest_n{ind};
        MSE=sum((yest-ydata).^2)/(length(ydata)-4); % 3 par + std err
        FPE=sum((yest-ydata).^2)*(length(ydata)+4)/(length(ydata)-4); % 3 par + std err
        SigmaSqr=std(yest-ydata)^2;
        LogLike=-0.5*length(ydata)*log(2*pi*SigmaSqr)-0.5*sum((yest-ydata).^2)/SigmaSqr;
        AIC=-2*LogLike+2*4;  % From Klaus Holst. 4 parameters is 3 model parameters and noise variance
        %ChiSqr=LogLike;
        %ChiSqr=sum((yest-ydata).^2)/std(ydata)^2; %From wikipedia
        %ChiSqr=sum((yest-ydata).^2)/std(yest-ydata)^2; %From wikipedia
        %ChiSqr=sum((yest-ydata).^2)/(std(yest-ydata).^2*(length(ydata)-3)); %From wikipedia
        %ChiSqr=sum((yest-ydata).^2)/(length(ydata)-3); %From wikipedia
        %ChiSqr=sum((yest-ydata).^2)/(7200/5*1000); %From wikipedia
        %AIC=ChiSqr+2*3;  %From wikipedia, 3 is no of parameters
    end


    function Ct=SRTM2(p,ref)
        %
        % Model of Reference tissue solution  (eq 6)
        %
        R1=p(1);
        BPnd=p(2);
        %
        % C(t)=R1 * C'(t) + R1 * [k2' - k2a] * C'(t) conv exp(-k2a*t)
        %
        %
        % BPnd=R1*k2'/k2a-1;
        %
        % k2a=R1*k2'/(BPnd+1)
        % b1=R1; b2=R1 *[k2' - k2a];
        %
        %
        k2a=R1*k2p/(BPnd+1);
        Ct=R1*ref(:,2)+R1*(k2p-k2a)*exp(-k2a*ref(:,1)).*KinmodCumtrapz_l(ref(:,1),ref(:,2).*exp(k2a*ref(:,1)));
        %
    end

    function [ix,t]=spline_int(t,x)
        %
        % Spline interpolation of function and quadgk integration of this
        %
        pp=spline(t,x);
        ix=quadgk(@(t) ppval(pp,t), min(t), max(t));
    end

end



