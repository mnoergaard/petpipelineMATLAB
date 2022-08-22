function ix=KinmodCumtrapz_l(t,x)
%
% Integration using rectangle rule with varying time steps
%
% It has been optimized for the varying time steps, requires that
% the time "t" is given as time in middle of sampling interval, and
% that sampling starts to time 0. The value x should also be the
% average value in the time interval.
%
% Like (start timing at time=0sec):
%   t (sec)    x
%   5          10
%   15         200
%   35         600
%
%
% Integration using trapezoidal rule with varying time steps
%
% CS, 20180423
% CS, 20200324
%

%
% Since 2018 a trapz method has been used instead of cumsum.
% Implemented by CS, to take care of first step with knowledge about times
% being in the middle of the sampling interval
% Also mathworks implementation using cumtrapz has been tried out
%
% Unfortunately for the ESRTM method (March 2020)it has been seen that for
% the first time after the step, the trapz method doesn't work properly,
% gives a overshoot in the estimate of the integral, therefore we have
% switched back to use the cumsum method
%

% CS implementation of trapz, does not work properly for ESRTM
% dt=DeltaT(t);
% tmp1=((x(2:end)-x(1:end-1))/2+x(1:end-1)).*(dt(1:end-1)+dt(2:end))/2;
% ix=cumsum([x(1)/2*dt(1)/2;tmp1]);


% Mathworks implementation of trapz, works as well as cumsum
% ix=cumtrapz(t,x);

% CS implemantation of cumsum (same as used before 2018)
dt=DeltaT(t);
tmp1=((x(2:end)-x(1:end-1))*1/2+x(1:end-1))*1/2.*(dt(1:end-1)+dt(2:end));
ix=cumsum([1/2*dt(1)*x(1);tmp1]);

    % Part of the CS implementation of trapz
%     function dt=DeltaT(t)
%         dt=zeros(size(t));
%         for k=1:length(t)
%             if k==1
%                 dt(k)=(t(2)-t(1));
%             else
%                 dt(k)=(t(k)-(t(k-1)+dt(k-1)/2))*2;
%             end
%         end
%     end

    % CS original cumsum implementation
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