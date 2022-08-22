function [Err,Par,yest,Par_cov]=frtm_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC,Name,Iter)
%
%  function  [Err,Par,yest]=frtm_estm(TimeTAC,FrameWeight,RefTAC,RoiTAC[,Iter])
%
% Implementation of FRTM estimation
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
%
%   Par       - [R1, k2, k3, BPnd]
%
% CS, 20140318
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
    [Err(i), Par{i}, yest(:,i), Par_cov{i}]=EstmRefTissue;
end
%

    function [Err, Par, yest, Par_cov]=EstmRefTissue
        %
        %
        for k=1:Iter
            fprintf('  Iteration: %i',k);
            %
            Par0=[1,0.5,0.5,2];   % R1, k2, k3, BPnd
            Par0=(1+0.1*randn(1,4)).*Par0;
            optims.MaxIter=1000;
            opt=statset;
            opt.FunValCheck='off';
            %[Par_n{k},R,J,Par_cov_n{k},Err_n(k)]=nlinfit(xdata(:,1:2),ydata,@RefTissue,Par0,opt);
            [Par_n{k},R,J,Par_cov_n{k},Err_n(k)]=nlinfit(xdata(:,1:2),ydata,@FRTM,Par0,opt,'Weights',xdata(:,3));
            yest_n{k}=FRTM(Par_n{k},xdata(:,1:2));
            fprintf(', LSQ: %6.3e\n',Err_n(k));
        end
        [minv,ind]=min(Err_n);
        Par=[Par_n{ind}];
        Par_cov=[Par_cov_n{ind}];
        Err=Err_n(ind);
        yest=yest_n{ind};
    end


    function Ct=FRTM(p,ref)
        %
        % Model of Reference tissue solution 
        Impl=1;
        %
        if Impl==1
            %
            % Implementeret som beskrevet i PET Pharmacokinetic Book
            %
            % Comvolution implementeret efter egen ide
            %
            R1=p(1);
            k2=p(2);
            k3=p(3);
            BP=p(4);
            %
            k4=k3/BP;
            %
            s=k2+k3+k4;
            r=k2/R1;
            q=4*k2*k4;
            p=sqrt(s^2-q);
            d=(s-p)/2;
            c=(s+p)/2;
            b=(d-k3-k4)*(d-r)/p;
            a=(k3+k4-c)*(c-r)/p;
            %
            TrzMth=1;
            %
            if TrzMth==1
                %
                %             for j=1:size(ref,1)
                %                 tmpsum1=cumtrapz_l(ref(1:j,1),ref(1:j,2).*exp(-c*(ref(j,1)-ref(1:j,1))));
                %                 sum1(j,1)=tmpsum1(end);
                %                 tmpsum2=cumtrapz_l(ref(1:j,1),ref(1:j,2).*exp(-d*(ref(j,1)-ref(1:j,1))));
                %                 sum2(j,1)=tmpsum2(end);
                %             end
                sum1=zeros(size(ref,1),1);
                sum2=zeros(size(ref,1),1);
                for j=1:size(ref,1)
                    sum1(j)=trapz_l(ref(1:j,1),ref(1:j,2).*exp(-c*(ref(j,1)-ref(1:j,1))));
                    sum2(j)=trapz_l(ref(1:j,1),ref(1:j,2).*exp(-d*(ref(j,1)-ref(1:j,1))));
                end
                %
                Ct=R1*(ref(:,2)+a*sum1+b*sum2);
                %
            else
                tmp=zeros(size(ref(:,1)));
                tmp(1)=1/2*1/2*ref(1,1)*ref(1,2).*exp(P0*ref(1,1));
                %
                % Fra hjemmesiden:
                % http://stackoverflow.com/questions/13396546/numerical-integration-over-non-uniform-grid-in-matlab-is-there-any-function
                %
                pp=spline(ref(:,1),ref(:,2).*exp(P0*ref(:,1)));
                %
                for k=2:length(ref(:,1))
                    %tmp(k)=tmp(1)+spline_int(ref(1:k,1),ref(1:k,2).*exp(P0*ref(1:k,1)));
                    tmp(k)=tmp(1)+quadgk(@(t) ppval(pp,t), ref(1,1), ref(k,1));
                end
                Ct=R1*ref(1:k,2)+(k2-R1*P0)*exp(-P0*ref(1:k,1)).*tmp;
            end
        else
            %
            % Implementeret som beskrevet i Lammertsma 1996 artikel
            %
            % Comvolution implementeret efter egen ide
            %
            K1=p(1);
            k2=p(2);
            k3=p(3);
            BP=p(4);
            %
            k4=k3/BP;
            %
            s=k2+k3+k4;
            q=4*k2*k4;
            p=sqrt(s^2-q);
            r=K1/p;
            d=(s+p)/2;
            c=(s-p)/2;
            b=(d-k3-k4)*r;
            a=(k3+k4-c)*r;
            %
            sum1=zeros(size(ref,1),1);
            sum2=zeros(size(ref,1),1);
            for j=1:size(ref,1)
                sum1(j)=trapz_l(ref(1:j,1),ref(1:j,2).*exp(-c*(ref(j,1)-ref(1:j,1))));
                sum2(j)=trapz_l(ref(1:j,1),ref(1:j,2).*exp(-d*(ref(j,1)-ref(1:j,1))));
            end
            %
            Ct=a*sum1+b*sum2;
            %
        end
    end

    function [ix,t]=spline_int(t,x)
        %
        % Spline interpolation of function and quadgk integration of this
        %
        pp=spline(t,x);
        ix=quadgk(@(t) ppval(pp,t), min(t), max(t));
    end


    function ix=cumtrapz_l(t,x)
        %
        % Integration using rectangle rule with varying time steps
        %
        %dt=DeltaT(t);
        %dx=dt.*x;
        %ix=cumsum(dx)-0.5*dx;
        %
        % Integration usingtrapezoidal rule with varying time steps
        %
        dt=DeltaT(t);
        tmp1=((x(2:end)-x(1:end-1))*1/2+x(1:end-1))*1/2.*(dt(1:end-1)+dt(2:end));
        tmp2=cumsum(tmp1);
        ix=[1/2*dt(1)*x(1);tmp2];
    end

    function ix=trapz_l(t,x)
        %
        % Integration using rectangle rule with varying time steps
        %
        %dt=DeltaT(t);
        %dx=dt.*x;
        %ix=cumsum(dx)-0.5*dx;
        %
        % Integration usingtrapezoidal rule with varying time steps
        %
        dt=DeltaT(t);
        tmp1=((x(2:end)-x(1:end-1))*1/2+x(1:end-1))*1/2.*(dt(1:end-1)+dt(2:end));
        tmp1=[1/2*dt(1)*x(1);tmp1];
        ix=sum(tmp1);
    end


    function dt=DeltaT(t)
        dt=zeros(size(t));
        for k=1:length(t)
            if k==1
                dt(k)=t(k)*2;
            else
                dt(k)=(t(k)-(t(k-1)+1/2*dt(k-1)))*2;
            end
        end
    end
end



