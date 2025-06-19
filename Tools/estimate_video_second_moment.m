function [sigObs,sigOerr] = estimate_video_second_moment(xObs,transect_video,SEtransect_video,mu0,xb) 
    dx=xObs(2)-xObs(1);
    [~,ixb]=min(abs(xObs-xb));
    xObs=xObs(ixb:end);

    % 1st moment
    muObs=[];
    for t=1:size(transect_video,2)
        tmp1=trapz(xObs,squeeze(transect_video(ixb:end,t)));
        muObs(t)=trapz(xObs,squeeze(transect_video(ixb:end,t)').*xObs)./tmp1;
    end
    if mu0~=99
        muObs=muObs*0+mu0;
    end
    
    % 2nd moment
    sigObs=[];
    sigOerr=[];
    for t=1:size(transect_video,2)
        SE=SEtransect_video(ixb:end,t);
        B=sum(transect_video(ixb:end,t)*dx);
        A=sum((xObs-muObs(t)).^2.*transect_video(ixb:end,t)'*dx);
        Delta_A = sum(xObs.^2 .* SE')*dx;
        Delta_B = sum(SE)*dx;
        sigObs(t)=A./B;
        sigOerr(t)=sqrt((Delta_A/B)^2 + (A*Delta_B/B^2)^2);
    end
end