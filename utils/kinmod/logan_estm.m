function [Err,Par,yest,Par_cov,LogLike,AIC, MSE, FPE]=logan_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name,Iter)
%
%  function  [Err,Par,yest]=logan_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name[,Iter])
%
% Implementation of Logan estimation
%
%   TimeTAC   - mid-time for each frame (in min, assumes that frame 1 starts at
%               time zero, so e.g. if length of first frame is 1/6 min, then
%               first time point should be 1/12 min)
%   FrameWeight - Weight of frame when fitting:
%                   W=(FrameLength^2/Trues)*exp_correct_factor
%   RefTAC    - Reference tissue TAC
%   RoiTAC    - Roi tissue TAC (can contain multiple tissue curves as columns)
%   Name      - Name.Ref - String with name of ref roi
%               Name.Roi - Cell array with names of regions
%               Name.Label - Label string for data series 
%   Iter      - No of fresh startup's (default is 10)
%
%   Par       - [R1, k2, BPnd]
%
% CS, 20150522
%
%
if nargin<6
    Iter=10;
end
%
xdata=[TimeTAC RefTAC FrameWeight];
%
for i=1:size(RoiTAC,2)
    ydata=RoiTAC(:,i);
    fprintf('\nEstimating ROI (%i/%i): %s\n',i,size(RoiTAC,2),Name.Roi{i});
    [Err(i), Par{i}, yest(:,i), Par_cov{i}, LogLike(i), AIC(i), MSE(i), FPE(i)]=EstmLoganTissue;
end
%

    function [Err, Par, yest, Par_cov, LogLike, AIC, MSE, FPE]=EstmLoganTissue
        %
        %
        % Original Logan formulation
        %
        % Int(Ct)./Ct = -1/Kappa2 + Vt * Int(Ca)./Ct
        % Kappa2 = k4 * k2 / (K3+k4)
        % Vt = K1/k2 * (k3+k4)/k4
        % 
        StSampl=7;
        %
        Cref_int=KinmodCumtrapz_l(xdata(:,1),RefTAC);
        Croi_int=KinmodCumtrapz_l(xdata(:,1),ydata);
        %
        x=Cref_int./ydata;
        y=Croi_int./ydata;
        %
        x=x(end-StSampl:end,1);
        y=y(end-StSampl:end,1);
        %
        glm_weight=xdata(end-StSampl:end,3);
        %
        %
        % Linear glm model
        %
        [b,dev,stats] = glmfit(x, y, 'normal', 'constant', 'on', 'weights',glm_weight);
        %[b,dev,stats] = glmfit(x, y, 'normal', 'constant', 'off');
        lin.yhat=[ones(length(x),1) x]*b;
        %lin.Err=sum((y-lin.yhat).^2.*xdata(:,3))/(length(y)-3);
        lin.Err=sum((y-lin.yhat).^2)/(length(y)-3);

        DEBUG=0;
        if DEBUG==1
            figure
            plot(x(:,1),y,'*',x(:,1),lin.yhat)
            Kappa2=-1/b(1)
            Vt=b(2)
            title(sprintf('Logan [%i-%i min], Kappa2=%5.3e, Vt=%5.3e',round(xdata(end-StSampl,1)),round(xdata(end,1)),Kappa2,Vt))
            xlabel('\int Ca / Ct')
            ylabel('\int Ct / Ct')
            print -dpsc2 -append Logan2impl.ps
        end
        
        lin.Kappa2=-1/b(1);
        lin.Vt=b(2);
        %lin.Kappa2_var=1/stats.covb(1,1);
        lin.Kappa2_var=NaN;
        lin.Vt_var=stats.covb(2,2);
        %
        Par=[lin.Kappa2; lin.Vt];
        Par_cov(1,1)=lin.Kappa2_var;
        Par_cov(2,2)=lin.Vt_var;
        Err=lin.Err;
        yest=lin.yhat;
        %
        % Calc stats
        %
        % Fast solution for error calculation, only using data points
        % included, CS, 20170208
        ydata=y;
        %
        MSE=sum((yest-ydata).^2)/(length(ydata)-3); % 2 par + std err
        FPE=sum((yest-ydata).^2)*(length(ydata)+3)/(length(ydata)-3); % 2 par + std err
        SigmaSqr=std(yest-ydata)^2;
        LogLike=-0.5*length(ydata)*log(2*pi*SigmaSqr)-0.5*sum((yest-ydata).^2)/SigmaSqr;
        AIC=-2*LogLike+2*4;  % From Klaus Holst. 4 parameters is 3 model parameters and noise variance
        %
        %pause
    end

end



