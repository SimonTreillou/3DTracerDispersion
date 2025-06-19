%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Figure 17 of Treillou et al. (2025)
% Code by Simon Treillou
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
addpath(genpath("Tools"));
addpath(genpath("/Users/simon/Code/CONFIGS/SURFZONE_MIXING"));
addpath(genpath("/Users/simon/Code/CONFIGS/IB09_SZ"));
%================== User defined parameters ===========================
%
% --- model params ---
%
% Simulations names and initial time of the release (s)
fnames{1}='IB09_3D_ShD'; tinits{1}=0;
fnames{2}='IB09_3D_MR_Dcrit'; tinits{2}=1800;
fnames{3}='IB09_3D_TRC_Dcrit_WMPER_test'; tinits{3}=1800;
fnames{4}='IB09_2D_BD'; tinits{4}=0;
fnames{5}='IB09_2D_TRC_Dcrit_WMPER_test'; tinits{5}=1800;
fnames{6}='IB09_2D_BD_TRC_Dcrit_WMPER_testDB15'; tinits{6}=1800;
fnames{7}='IB09_2D_BD3_TRC_Dcrit_WMPER_testDB15'; tinits{7}=1800;

% Simulation to display
display = [1,2,3,4,5,6,7];
[c1,c2]=define_colors(0.95);
% integ: depth-integrated concentration - precise: concentration at defined
% level
whlevel="integ"; 
makepdf=true;
names = ["3D - Shear Dis.", "3D - Mini-rips", "3D - Flash rips", ...
         "2D - $\kappa_0^a$", "2D - Flash rips ($\kappa_0$=0)", "2D - Flash rips ($\kappa_0$=0.075)", "2D - Flash rips ($\kappa_0$=0.3)"];
% Ensemble=True -> uses ensemble-averaged results for 3D-TRC
ensemble=true;
% Diffusion: "IS+SZ" for whole domain diffusivity; "SZ" for SZ diffusivity
diff="IS+SZ";
%======================================================================

%% MODEL 
% Init
sigs = cell(1,length(display));
sigerrs = cell(1,length(display));
mus  = cell(1,length(display));
tps  = cell(1,length(display));
cs   = cell(1,length(display));
ss   = cell(1,length(display));
kxxs = cell(1,length(display));
kxxstd  = cell(1,length(display));
n=1;
for i=display
    fname=[fnames{i},'/rip_avg.nc'];
    if i==1
        c=c1;
        s=':';
    elseif i==2
        c=c1;
        s="--";
    elseif i==3
        c=c1;
        s="-";
    elseif i==4
        c=[0.5,0.5,0.5];
        s="-";
    elseif i==5
        c=c2;
        c(4)=0.2;
        s="--";
    elseif i==6
        c=c2;
        c(4)=0.5;
        s="-.";
    elseif i==7
        c=c2;
        s="-";
    elseif i==8
        c=[0.3,0.3,0.9];
    elseif i==9
        c=[0.3,0.7,0.9];
    elseif i==10
        c=[0.3,0.7,0.2];
    end

    % Load files
    nc=netcdf(fname);

    % Time
    time=nc{'scrum_time'}(:);
    [~,tstr]=min(abs(time-tinits{i}));
    [~,tend]=min(abs(time-time(tstr)-6000));
    time=time(tstr:tend)-time(tstr);
    tp=time-time(1);

    % Bathymetry
    yr=squeeze(nc{'y_rho'}(:,1));
    xr=findx_IB09(nc);
    dx=xr(2)-xr(1);

    % Set cross-shore length
    if diff=="IS+SZ"
    % Lx is the cross-shore domain length. When considering mini-rips,
    % reduced to 120 m (~99% of cross-shore mass extent) to stay in the
    % linear diffusion regime
        Lx=300;
        if i==2
            Lx=120; 
        end
    elseif diff=="SZ"
        Lx=80; %% surfzone cross-shore length
    end
    [~,xb] = min(abs(xr+Lx));
    xmax=0;
    [~,x0] = min(abs(xr+xmax));
    xr=xr(xb:x0);

    % Tracer concentration depending on the desired level
    %  - integ: depth integration
    %  - precise: surface
    if whlevel=="integ"
        hr=squeeze(nc{'h'}(1,xb:x0));
        Dcrit=nc{'Dcrit'}(:);
        N=length(nc('s_rho'));
        theta_s=nc.theta_s(:); 
        theta_b=nc.theta_b(:); 
        hc=nc.hc(:); 
        zeta=squeeze(nc{'zeta'}(end,1,xb:x0));
        zw=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'w',2));
        dz=zw(2:end,:)-zw(1:end-1,:);         % ---> zw(2:N,:)
        t=squeeze(nc{'tpas01'}(tstr:tend,:,:,xb:x0));
        t=squeeze(sum(dz.*permute(t,[2 4 1 3]))./mean(hr));
        t=permute(t,[2 3 1]);

        if i==3
            nc2=netcdf([fnames{i},'2/rip_avg.nc']);
            zeta=squeeze(nc2{'zeta'}(end,1,xb:x0));
            zw=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'w',2));
            dz=zw(2:end,:)-zw(1:end-1,:);         % ---> zw(2:N,:)
            t2=squeeze(nc2{'tpas01'}(:,:,:,xb:x0));
            t2=squeeze(sum(dz.*permute(t2,[2 4 1 3]))./mean(hr));
            t2=permute(t2,[2 3 1]);
            t=cat(1,t,t2);
        end
    elseif whlevel=="precise"
         t=squeeze(nc{'tpas01'}(tstr:tend,10,:,xb:x0));
    end
    
    if i==3
        time=nc{'scrum_time'}(:);
        time=cat(1,time,nc2{'scrum_time'}(:)); 
        [~,tstr]=min(abs(time-tinits{i}));
        [~,tend]=min(abs(time-time(tstr)-6000));
        time=time(tstr:tend)-time(tstr);
        tp=time-time(1);
    end
    % Computation of 1st moment
    mu=[];
    for tt=1:(tend-tstr+1)
        mu(tt)=C10_mu(squeeze(t(tt,:,:)),xr);
    end
    mus{i}=mu;
    mu=mu*0; %absolute diffusion

    sig=[];sigerr=[];
    for tt=1:(tend-tstr+1)
        Dbar=squeeze(mean(t(tt,:,:),2));
        SE=squeeze(std(t(tt,:,:))/sqrt(size(t,2)));

        B=sum(Dbar*dx);
        A=sum(((xr-mu(tt)).^2.*Dbar')');
        
        Delta_A = sum((xr-mu(tt)).^2 .* SE');
        Delta_B = sum(SE);


        sig(tt)=A/B;
        sigerr(tt)= sqrt((Delta_A/B)^2 + (A*Delta_B/B^2)^2);
    end
    if i==3 && ensemble
        sig=importdata("./Data/Ensemble/3D_ensemble_sig_mean"+string(Lx)+".txt")';
        %sig=sig(1:length(time));
        sigstd=importdata("./Data/Ensemble/3D_ensemble_sig_stderr"+string(Lx)+".txt")';
        %sigstd=sigstd(1:length(time));
        time=time(1:length(sig));  
    end
    sigs{i}=sig-min(sig);
    sigerrs{i}=sigerr;

    % Computation of saturated 2nd moment
    sat=((xr(end)^3-xr(1)^3)/3) / (xr(end)-xr(1));
    R=sig/sat;
    R0=0.55;
    R=(R(2:end)+R(1:end-1))/2;
    if i==1
        R1=R;
    elseif i==2
        R2=R;
    end
    for j=1:length(R)
        if R(j)>R0
            break
        end
    end
    if j>229
        j=229;
    end
    ixR0=1:j;
    disp(j);

    % Computation of diffusivity
    dt=time(2)-time(1);
    kxx=0.5*((sig(2:end)-sig(1:end-1)))'./dt;
    n=length(sig(ixR0));
    a = sum((time(ixR0)-mean(time(ixR0))).*(sig(ixR0)'-mean(sig(ixR0)))) ./ sum((time(ixR0)-mean(time(ixR0))).^2);
    b = mean(sig(ixR0)')-a*mean(time(ixR0));
    errkxx = sqrt(mean(sigerr.^2) / sum((time(ixR0)-mean(time(ixR0))).^2));
    r2=1-sum((sig(ixR0)'-time(ixR0)*a-b).^2)./sum((sig(ixR0)'-mean(sig(ixR0))).^2);

    kxxs{i}=a*0.5;
    kxxstd{i}=errkxx*0.5;
    tps{i}=time-time(1);
    cs{i}=c;
    ss{i}=s;
    ixR0s{i}=ixR0;
    disp("------------------------------------")
    disp(names(i)+" -----> "+num2str(kxxs{i}));
    disp(r2);
    disp(errkxx);
    n=n+1;
end

%% PLOT
figure('Position',[100 100 900 500],'PaperOrientation','landscape');

displayplot=[1,2,3,4,6,7];
legs=cell(1,length(displayplot));
for i=1:length(displayplot)
    j=displayplot(i);
    tp=tps{j};
    ixR0=ixR0s{j};
    esp=2.5;
    if i==3
        esp=3.5;
    end
    plot(tp,smooth(sigs{j},5),'Color',cs{j},'LineWidth',esp,'LineStyle',ss{j});
    hold on
    scatter(tp(ixR0(end)),sigs{j}(ixR0(end)),1000,cs{j}(1:3),'.')
    if i==3 && ensemble
        patch([tp',fliplr(tp')], [smooth(sigs{j}-sigstd,5)',fliplr(smooth(sigs{j}+sigstd,5)')], cs{j}(1:end-1),'FaceAlpha',.3,'EdgeColor','none');
    end
    C=cs{j};C(end)=0.3;
    legs{i}=names{j}+" ($\kappa_{xx} = $ "+string(round(kxxs{j},2))+")";
end
legs={legs{1},'',legs{2},'',legs{3},'','',legs{4},'',legs{5},'',legs{6}};
ylabel("$\sigma^2$ ($m^2$)",'Interpreter','latex','FontSize',15);
xlabel("Time (s)",'Interpreter','latex','FontSize',15);
grid("minor");
legend([legs],'Interpreter','latex','FontSize',17,'Location','best');
xlim([0 5950]);
ylim([-10 15500])
set(gca,'linewidth',1.5);
set(gcf,'Color','w');


if makepdf
    print(gcf, '-dpdf', '-r600', ['./Figures/Figure17.pdf']);
end