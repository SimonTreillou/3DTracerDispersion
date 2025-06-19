function [rms] = RMSE_differentarrays(x1,y1,x2,y2)
    rms=0;
    if length(x1)<length(x2)
        for i=1:length(x1)
            [~,ix]=min(abs(x2-x1(i)));
            rms=rms+(y1(i)-y2(ix))^2;
        end
        rms=sqrt(rms/length(x1));
    elseif length(x1)>=length(x2)
        for i=1:length(x2)
            [~,ix]=min(abs(x1-x2(i)));
            rms=rms+(y2(i)-y1(ix))^2;
        end
        rms=sqrt(rms/length(x2));
    end
end