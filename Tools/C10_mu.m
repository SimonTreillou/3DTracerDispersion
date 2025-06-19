function [mu] = C10_mu(t,xr)
    % See Clark et al. 2010, eq. 6
    %INPUT:
    %  - t:  tracer concentration (time,y,x)
    %  - xr: cross-shore coordinates
    %OUTPUT:
    %  - mu: 1st moment
    
    Dbar=squeeze(mean(t));
    tmp1=trapz(xr',(Dbar)');
    mu=trapz(xr',(Dbar.*xr)')./tmp1;
end