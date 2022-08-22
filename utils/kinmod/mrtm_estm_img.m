function [Err,Par,yest,Par_cov,LogLike,AIC, MSE, FPE]=mrtm_estm_img(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name)
%
%  function  [Err,Par,yest]=mrtm_estm_img(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name)
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
%
%   Par -      [R1, k2', BPnd]   [k2=k2'*R1/(BPnd+1)]
%
% CS, 20190322
%
%
xdata=[TimeTAC RefTAC FrameWeight];
%
Par=zeros(size(RoiTAC,2),3);
Par_cov=zeros(size(RoiTAC,2),3,3);
for i=1:size(RoiTAC,2)
    ydata=RoiTAC(:,i);
    if any(isnan(ydata))
        fprintf('\nNot estimating ROI (%i/%i): %s, nan in data\n',i,size(RoiTAC,2),Name.Roi{i});
        Err(i)=nan;
        Par(i,:)=nan(1,3);
        yest(:,i)=nan(size(ydata));
        Par_cov(i,:,:)=nan(3,3);
        LogLike(i)=nan;
        AIC(i)=nan;
        MSE(i)=nan;
        FPE(i)=nan;
    else
        %fprintf('\nEstimating ROI (%i/%i): %s\n',i,size(RoiTAC,2),Name.Roi{i});
        [Err(i), Par(i,:), yest(:,i), Par_cov(i,:,:), LogLike(i), AIC(i), MSE(i), FPE(i)]=EstmMRTMTissue;
        %if (Par(i,3)<-1)&&(Par(i,3)>-2)
            %figure;plot(TimeTAC,RefTAC,TimeTAC,yest(:,i));
            %keyboard
        %end
    end
end
%

    function [Err, Par, yest, Par_cov, LogLike, AIC, MSE, FPE]=EstmMRTMTissue
        %
        %
        % From kinetic course notes
        %
        % Ct(T)=R1K2' * Int[0;T] C'(t)dt - k2 * Int[0;T] C(t)dt + R1 C'(T)
        % BPnd=-(b1/b2+1); R1=b3; k2'=b1/b3;
        %
        Cref_int=KinmodCumtrapz_l(xdata(:,1),RefTAC);
        Croi_int=KinmodCumtrapz_l(xdata(:,1),ydata);
        %
        x=[Cref_int,Croi_int,RefTAC];
        y=ydata;
        %
        %
        % Linear glm model
        %
        [b,dev,stats] = glmfit(x, y, 'normal', 'constant', 'off', 'weights',xdata(:,3));
        lin.yhat=x*b;
        lin.Err=sum((y-lin.yhat).^2.*xdata(:,3))/(length(y)-3);
        lin.R1=b(3);
        lin.k2p=b(1)/b(3);
        lin.BPnd=-(b(1)/b(2)+1);
        lin.R1_var=stats.covb(3,3);
        lin.k2p_var=(stats.covb(1,1)*b(3)^2-2*stats.covb(1,3)*b(1)*b(3)+stats.covb(3,3)*b(1)^2)/...
            b(3)^4;
        lin.BPnd_var=(stats.covb(1,1)*b(2)^2-2*stats.covb(1,2)*b(1)*b(2)+stats.covb(2,2)*b(1)^2)/...
            b(2)^4;
        Par=[lin.R1 lin.k2p lin.BPnd];
        Par_cov=zeros(3);
        Par_cov(1,1)=lin.R1_var;
        Par_cov(2,2)=lin.k2p_var;
        Par_cov(3,3)=lin.BPnd_var;
        Err=lin.Err;
        yest=lin.yhat;
        %
        % Calc stats
        %
        MSE=sum((yest-ydata).^2)/(length(ydata)-4); % 3 par + std err
        FPE=sum((yest-ydata).^2)*(length(ydata)+4)/(length(ydata)-4); % 3 par + std err
        SigmaSqr=std(yest-ydata)^2;
        LogLike=-0.5*length(ydata)*log(2*pi*SigmaSqr)-0.5*sum((yest-ydata).^2)/SigmaSqr;
        AIC=-2*LogLike+2*4;  % From Klaus Holst. 4 parameters is 3 model parameters and noise variance
        %
    end


    function Ct=MRTM(p,xd)
        %
        % Model of Reference tissue solution  (eq 6)
        %
        R1=p(1);
        k2p=p(2);
        BPnd=p(3);
        %
        % Ct(T)=R1K2' * Int[0;T] C'(t)dt - k2 * Int[0;T] C(t)dt + R1 C'(T)
        % BPnd=-(b1/b2+1); R1=b3; k2'=b1/b3;
        %
        % b3=R1; b1=k2'*R1; b2=-b1/(1+BPnd)=-k2'*R1/(1+BPnd)
        %
        Ct=k2p*R1*xd(:,1) - k2p*R1/(1+BPnd)*xd(:,2) + R1*xd(:,3);
        %

    end        
        
end



