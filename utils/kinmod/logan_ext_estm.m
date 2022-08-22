function [Err,Par,yest,Par_cov,LogLike,AIC, MSE, FPE]=logan_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name,TimeChlng)
%
%  function  [Err,Par,yest]=logan_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name)
%
% Implementation of Logan estimation, with a challenge, so last 7 point
% before challenge is used for estimation of slope before and last 7 point
% in time series used for estimation after challenge
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
%
%   Par       - [R1, k2, BPnd]
%
% CS, 20150522
%
%
%
xdata=[TimeTAC RefTAC FrameWeight];
%
%
% Calculate each frame finish time and then calculate challenge index
% (minimun time difference to challenge time)
%
for i=1:length(TimeTAC)
    if i==1
        FrameLength(i)=TimeTAC(1)*2;
    else
        FrameLength(i)=(TimeTAC(i)-TimeTAC(i-1)-FrameLength(i-1)/2)*2;
    end
    FrameFinishTime(i)=TimeTAC(i)+FrameLength(i)/2;
end
[TimeDiff,CH_ind]=min(abs(FrameFinishTime-TimeChlng));
fprintf('Challenge implemented at: %5.1f [min] (distance between frame finish and challenge time is %5.1f)\n',TimeTAC(CH_ind),TimeDiff);
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
        % Challenge model before and after
        %
        dd=CH_ind;
        %
        for k=0:1
            %
            if k==0 % Before challenge
                xfit=x(CH_ind-StSampl:CH_ind,1);
                yfit=y(CH_ind-StSampl:CH_ind,1);
                %
                glm_weight=xdata(CH_ind-StSampl:CH_ind,3);
                %
                Time.Start=xdata(CH_ind-StSampl,1);
                Time.Stop=xdata(CH_ind,1);
            else
                xfit=x(end-StSampl:end,1);
                yfit=y(end-StSampl:end,1);
                %
                glm_weight=xdata(end-StSampl:end,3);
                %
                Time.Start=xdata(end-StSampl,1);
                Time.Stop=xdata(end,1);
            end
            %
            % Linear glm model
            %
            [b,dev,stats] = glmfit(xfit, yfit, 'normal', 'constant', 'on', 'weights',glm_weight);
            %[b,dev,stats] = glmfit(x, y, 'normal', 'constant', 'off');
            
            if k==0
                lin.yhat_0=[ones(length(xfit),1) xfit]*b;
                lin.Err_0=sum((yfit-lin.yhat_0).^2)/(length(yfit)-3);
                %
                lin.Kappa2_0=-1/b(1);
                lin.Vt_0=b(2);
                %lin.Kappa2_var=1/stats.covb(1,1);
                lin.Kappa2_var_0=NaN;
                lin.Vt_var_0=stats.covb(2,2);
            else
                lin.yhat_1=[ones(length(xfit),1) xfit]*b;
                lin.Err_1=sum((yfit-lin.yhat_1).^2)/(length(yfit)-3);
                %
                lin.Kappa2_1=-1/b(1);
                lin.Vt_1=b(2);
                %lin.Kappa2_var=1/stats.covb(1,1);
                lin.Kappa2_var_1=NaN;
                lin.Vt_var_1=stats.covb(2,2);
            end
            
            DEBUG=0;
            if DEBUG==1
                figure
                if k==0
                    plot(xfit(:,1),yfit,'*',xfit(:,1),lin.yhat_0);
                else
                    plot(xfit(:,1),yfit,'*',xfit(:,1),lin.yhat_1);
                end
                Kappa2=-1/b(1)
                Vt=b(2)
                title(sprintf('Logan [%i-%i min], Kappa2=%5.3e, Vt=%5.3e',round(Time.Start),round(Time.Stop),Kappa2,Vt))
                xlabel('\int Ca / Ct')
                ylabel('\int Ct / Ct')
                print -dpsc2 -append Logan2impl.ps
            end
            
        end
        
        %
        Par=[lin.Kappa2_0; lin.Vt_0;lin.Kappa2_1; lin.Vt_1];
        Par_cov(1,1)=lin.Kappa2_var_0;
        Par_cov(2,2)=lin.Vt_var_0;
        Par_cov(3,3)=lin.Kappa2_var_1;
        Par_cov(4,4)=lin.Vt_var_1;
        Err=lin.Err_1;
        yest=lin.yhat_1;
        %
        % Calc stats
        %
        % Fast solution for error calculation, only using data points
        % included, CS, 20170208
        ydata=yfit;
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



