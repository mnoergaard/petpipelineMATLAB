function [Err,Par,yest,Par_cov,LogLike,AIC, MSE, FPE]=logan_ref_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name,ExtPar)
%
%  function
%  [Err,Par,yest]=logan_ref_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name,ExtPar)
%
% Implementation of Logan reference estimation
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
%   ExtPar    - ExtPar.k2ref - The value of k2' (for the reference region)
%                              [1/min]
%               ExtPar.tstar - from this datapoint and on used for
%                              estimating linear slope [min]
%
%   Par       - [1/kappa2, BPnd]  
%                  BPnd = Vt/Vt' -1
%                  kappa2 = k2 * k4 / (k3+k4)
%
% CS, 20180509
%
%
xdata=[TimeTAC RefTAC FrameWeight];
%
for i=1:size(RoiTAC,2)
    ydata=RoiTAC(:,i);
    fprintf('\nEstimating ROI (%i/%i): %s\n',i,size(RoiTAC,2),Name.Roi{i});
    [Err(i), Par{i}, yest(:,i), Par_cov{i}, LogLike(i), AIC(i), MSE(i), FPE(i)]=EstmLoganRefTissue;
end
%

    function [Err, Par, yest, Par_cov, LogLike, AIC, MSE, FPE]=EstmLoganRefTissue
        %
        %
        % Original Logan formulation
        %
        % Int(Ct)./Ct = -1/Kappa2 + Vt * Int(Ca)./Ct
        % Kappa2 = k4 * k2 / (K3+k4)
        % Vt = K1/k2 * (k3+k4)/k4
        %
        % Original Logan Ref formulation
        %
        % Int(Ct)./Ct = -1/Kappa2 + Vt/Vt' * (Int(Cref) + Cref/k2')./Ct
        % Kappa2 = k4 * k2 / (K3+k4)
        % BPnd = Vt/Vt' - 1
        % 
        index=find(xdata(:,1)>ExtPar.tstar);
        if (isempty(index))
            StartSampl=1;
        else
            StartSampl=min(index);
        end
        %
        Cref_int=KinmodCumtrapz_l(xdata(:,1),RefTAC);
        Croi_int=KinmodCumtrapz_l(xdata(:,1),ydata);
        %
        x=(Cref_int+RefTAC/ExtPar.k2ref)./ydata;
        y=Croi_int./ydata;
        %
        xfit=x(StartSampl:end,1);
        yfit=y(StartSampl:end,1);
        %
        glm_weight_fit=xdata(StartSampl:end,3);
        %
        %
        % Linear glm model
        %
        [b,dev,stats] = glmfit(xfit, yfit, 'normal', 'constant', 'on', 'weights',glm_weight_fit);
        %[b,dev,stats] = glmfit(x, y, 'normal', 'constant', 'off');
        lin.yhat=[ones(length(xfit),1) xfit]*b;
        %lin.Err=sum((y-lin.yhat).^2.*xdata(:,3))/(length(y)-3);
        lin.Err=sum((yfit-lin.yhat).^2)/(length(yfit)-3);

        DEBUG=1;
        if DEBUG==1
            figure
            plot(x,y,'*',xfit(:,1),yfit,'o',xfit(:,1),lin.yhat)
            Kappa2=-1/b(1);
            BPnd=b(2)-1;
            title(sprintf('Logan [%i-%i min], Kappa2=%5.3e, BPnd=%5.3e',round(xdata(StartSampl,1)),round(xdata(end,1)),Kappa2,BPnd))
            ylabel('\int Ct / Ct')
            xlabel('(\int Cref + Cref/k2'') / Ct')
            print -dpsc2 -append Logan2impl.ps
        end
        
        lin.Kappa2=-1/b(1);
        lin.BPnd=b(2)-1;
        %
        lin.Kappa2_var=stats.covb(1,1)/b(1)^4;
        lin.BPnd_var=stats.covb(2,2);
        %
        Par=[lin.Kappa2; lin.BPnd];
        Par_cov(1,1)=lin.Kappa2_var;
        Par_cov(2,2)=lin.BPnd_var;
        Err=lin.Err;
        yest=lin.yhat;
        %
        % Calc stats
        %
        % Fast solution for error calculation, only using data points
        % included, CS, 20170208
        ydata=y;
        %
        MSE=sum((yest-yfit).^2)/(length(yfit)-3); % 2 par + std err
        FPE=sum((yest-yfit).^2)*(length(yfit)+3)/(length(yfit)-3); % 2 par + std err
        SigmaSqr=std(yest-yfit)^2;
        LogLike=-0.5*length(yfit)*log(2*pi*SigmaSqr)-0.5*sum((yest-yfit).^2)/SigmaSqr;
        AIC=-2*LogLike+2*4;  % From Klaus Holst. 4 parameters is 3 model parameters and noise variance
        %
        %pause
    end

end



