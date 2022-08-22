function [Err,Par,yest,Par_cov,LogLike,AIC,MSE,FPE]=esrtm_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name,TimeChlng,Iter)
%
%  function  [Err,Par,yest,Par_cov,ChiSqr,AIC]=esrtm_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name,TimeChlng[,Iter])
%
% Implementation of ESRTM estimation
%
%   TimeTAC   - mean time for each frame (in min, assumes that frame 1 starts at
%               time zero, so e.g. if length of first frame is 1/6 min, then
%               first time point should be 1/12 min)
%   FrameWeight - Weight of frame when fitting:
%                   W=(FrameLength^2/Trues)*exp_correct_factor
%   RefTAC    - Reference tissue TAC
%   RoiTAC    - Roi tissue TAC (can contain multiple tissue curves as columns)
%   Name      - Name.Ref - String with name of ref roi
%               Name.Roi - Ceaal array with names of regions
%               Name.Label - Label string for data series
%   TimeChlng - Time for challenge, where shift between two models is
%               performed
%
%   Par       - [R1, k2, BPnd0, BPnd1]
%
% CS, 20140120
%
%
if nargin<7
    Iter=10;
end
%
xdata=[TimeTAC RefTAC FrameWeight];
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
    if any(isnan(ydata))
        fprintf('\nNot estimating ROI (%i/%i): %s, nan in data\n',i,size(RoiTAC,2),Name.Roi{i});
        Err(i)=nan;
        Par{i}=nan(4,1);
        yest(:,i)=nan(size(ydata));
        Par_cov{i}=nan(4,4);
        LogLike(i)=nan;
        AIC(i)=nan;
        MSE(i)=nan;
        FPE(i)=nan;
    else
        fprintf('\nEstimating ROI (%i/%i): %s\n',i,size(RoiTAC,2),Name.Roi{i});
        [Err(i), Par{i}, yest(:,i), Par_cov{i}, LogLike(i), AIC(i), MSE(i), FPE(i)]=EstmExtRefTissue;
    end
end
%

    function [Err, Par, yest, Par_cov, LogLike, AIC, MSE, FPE]=EstmExtRefTissue
        %
        %
        for k=1:Iter
            fprintf('  Iteration: %i',k);
            Par0=(1+0.1*randn(1,4)).*[1,0.1,2,2];   % R1, k2, BP0, BP1
            optims.MaxIter=1000;
            opt=statset;
            opt.FunValCheck='off';
            [Par_n{k},R,J,Par_cov_n{k},Err_n(k)]=nlinfit(xdata(:,1:2),ydata,@ExtRefTissue,Par0,opt,'Weights',xdata(:,3));
            yest_n{k}=ExtRefTissue(Par_n{k},xdata(:,1:2));
            fprintf(', LSQ: %6.3e\n',Err_n(k));
        end
        [minv,ind]=min(Err_n);
        Par=[Par_n{ind}];
        Par_cov=[Par_cov_n{ind}];
        Err=Err_n(ind);
        yest=yest_n{ind};
        MSE=sum((yest-ydata).^2)/(length(ydata)-4); % 4 par + std err
        FPE=sum((yest-ydata).^2)*(length(ydata)+5)/(length(ydata)-5); % 4 par + std err
        SigmaSqr=std(yest-ydata)^2;
        LogLike=-0.5*length(ydata)*log(2*pi*SigmaSqr)-0.5*sum((yest-ydata).^2)/SigmaSqr;
        AIC=-2*LogLike+2*5;  % From Klaus Holst. 5 parameters is 4 model parameters and noise variance 
        %ChiSqr=LogLike;
%        ChiSqr=sum((yest-ydata).^2)/std(ydata)^2; %From wikipedia
%        AIC=ChiSqr+2*4;  %From wikipedia, 4 is no of parameters in ESRTM
    end

    function [C,yest]=CostExtRefTissue(Parm)
        %
        % Returns value of cost func
        %
        yest=ExtRefTissue(Parm,xdata);
        %
        %
        %C=sqrt(sum((ydata-yest).^2)/length(yest));
        C=sqrt(sum((ydata-yest).^2.*xdata(:,3))/length(yest));        
    end

    function Ct=ExtRefTissue(p,ref)
        %
        % Model of Reference tissue solution  (eq 6)
        %
        R1=p(1);
        k2=p(2);
        BP0=p(3);
        BP1=p(4);
        %
        P0=k2/(1+BP0);
        %
        Ct=zeros(size(ref(:,1)));
        %
        %
        TrzMth=1;
        %
        if TrzMth==1
            %
            Ct(1:CH_ind)=R1*ref(1:CH_ind,2)+(k2-R1*P0)*exp(-P0*ref(1:CH_ind,1)).*...
                KinmodCumtrapz_l(ref(1:CH_ind,1),ref(1:CH_ind,2).*exp(P0*ref(1:CH_ind,1)));
            %
            T0=ref(CH_ind,1);
            Ct_T0=Ct(CH_ind);
            Cref_T0=ref(CH_ind,2);
            P1=k2/(1+BP1);
            %
            Ct(CH_ind:size(ref,1))=Ct_T0*exp(-P1*(ref(CH_ind:end,1)-T0))+...
                R1*(ref(CH_ind:end,2)-Cref_T0*exp(-P1*(ref(CH_ind:end,1)-T0)))+...
                (k2-R1*P1)*exp(-P1*ref(CH_ind:end,1)).*...
                KinmodCumtrapz_l(ref(CH_ind:end,1)-T0,ref(CH_ind:end,2).*exp(P1*ref(CH_ind:end,1)));
            %
            Ct=reshape(Ct,size(ref(:,1)));
            %
        else
            tmp=zeros(size(ref(1:CH_ind,1)));
            tmp(1)=1/2*1/2*ref(1,1)*ref(1,2).*exp(P0*ref(1,1));
            %
            pp=spline(ref(:,1),ref(:,2).*exp(P0*ref(:,1)));
            %
            for k=2:length(ref(1:CH_ind,1))
                %tmp(k)=tmp(1)+spline_int(ref(1:k,1),ref(1:k,2).*exp(P0*ref(1:k,1)));
                tmp(k)=tmp(1)+quadgk(@(t) ppval(pp,t), ref(1,1), ref(k,1));
            end
            %
            Ct(1:CH_ind)=R1*ref(1:CH_ind,2)+(k2-R1*P0)*exp(-P0*ref(1:CH_ind,1)).*...
                tmp;
            %
            T0=ref(CH_ind,1);
            Ct_T0=Ct(CH_ind);
            Cref_T0=ref(CH_ind,2);
            P1=k2/(1+BP1);
            %
            tmp=zeros(size(ref(CH_ind:end,1)));
            tmp(1)=0;
            %
            pp=spline(ref(:,1),ref(:,2).*exp(P1*ref(:,1)));
            %
            for k=2:length(ref(CH_ind:end,1))
                %tmp(k)=spline_int(ref(CH_ind:CH_ind+k-1,1)-T0,ref(CH_ind:CH_ind+k-1,2).*exp(P1*ref(CH_ind:CH_ind+k-1,1)));
                tmp(k)=tmp(1)+quadgk(@(t) ppval(pp,t), ref(CH_ind,1), ref(CH_ind+k-1,1));
            end
            %
            Ct(CH_ind:size(ref,1))=Ct_T0*exp(-P1*(ref(CH_ind:end,1)-T0))+...
                R1*(ref(CH_ind:end,2)-Cref_T0*exp(-P1*(ref(CH_ind:end,1)-T0)))+...
                (k2-R1*P1)*exp(-P1*ref(CH_ind:end,1)).*...
                tmp;
        end
    end


    function [ix,t]=spline_int(t,x)
        %
        % Spline interpolation of function and quadgk integration of this
        %
        pp=spline(t,x);
        ix=quadgk(@(t) ppval(pp,t), min(t), max(t));
    end

end



