function [x] = findx_IB09(nc)
    xl=nc{'xl'}(1);
    h=nc{'h'}(1,:);
    [~,ix]=min(abs(h));
    xr=squeeze(nc{'x_rho'}(1,:));
    x=xr-xl+(xl-xr(1,ix));
end