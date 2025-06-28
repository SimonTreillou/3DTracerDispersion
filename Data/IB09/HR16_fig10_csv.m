%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%close all
%================== User defined parameters ===========================
%
% --- model params ---
%
%fname     = '../rip_his.nc';  % CROCO file name
%dirpath    = '/Users/simon/Code/IB09/IB09_PSbienpos/';
%dirpath    = '/Users/simon/Code/IB09/IB09PATdxbathynewang11strat10radiatif_newWM2023/';
dirpath    = '/Users/simon/Code/IB09/IB09_PSbienposT2/';
dirpath    = '/Users/simon/Code/CONFIGS/IB09_randomphase_S30/';
dirpath    = '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/';   % croco history file name
%dirpath    = '/Users/simon/Code/CONFIGS/IB09_S102/';   % croco history file name
%dirpath    = '/Users/simon/Code/Matlab/';
%dirpath    = '/Users/simon/Code/CONFIGS/IB09/';
%dirpath    = '/Users/simon/Code/CONFIGS/Compile_Patrick/'
fname      = strcat(dirpath,'rip_avgtot.nc');
%fname     = '../../../IB09_exp_Qnew_wind/rip_his.nc';

model3D   = 1;

makepdf   = 0;             % make pdf file
Dcrit=0.2;
%
%======================================================================

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);

dx=3;

tstr=10;

dyn="3D";
if dyn=="3D"
    load /Users/simon/Code/CONFIGS/IB09_S10_newPS/IB09_S10_newPS_time.csv
    time=IB09_S10_newPS_time;
    load /Users/simon/Code/CONFIGS/IB09_S10_newPS/IB09_S10_newPS_2_time.csv
    time=cat(1,time,IB09_S10_newPS_2_time);
    load /Users/simon/Code/CONFIGS/IB09_S10_newPS/IB09_S10_newPS_3_time.csv
    time=cat(1,time,IB09_S10_newPS_3_time);
elseif dyn=="2D"
    load /Users/simon/Code/CONFIGS/IB09_DB1tsdiff/IB09_DB1tsdiff_time.csv
    time=IB09_DB1tsdiff_time;
    load /Users/simon/Code/CONFIGS/IB09_DB1tsdiff/IB09_DB1tsdiff_2_time.csv
    time=cat(1,time,IB09_DB1tsdiff_2_time);
end
time=time(1:end-4);
tend=length(time);
if tend<tstr; tstr=1; end;
net_v=2.13e6;

h=nc{'h'}(1,:);
beach_l = 50;%abs(min(h)/0.02); % 0.02 is beach slope (cf HR16)
xl=nc{'xl'}(:);
x=nc{'x_rho'}(:,:)-xl+beach_l;
y=nc{'y_rho'}(:,:);
pm=nc{'pm'}(:,:);
pn=nc{'pn'}(:,:);
N=length(nc('s_rho'));

if model3D,
  ui=squeeze(nc{'u'}(1,N,:,:));
  vi=squeeze(nc{'v'}(1,N,:,:));
  u=squeeze(nc{'u'}(tstr:tend,N,:,:));
  v=squeeze(nc{'v'}(tstr:tend,N,:,:));
  ub=squeeze(nc{'u'}(tstr:tend,1,:,:));
  vb=squeeze(nc{'v'}(tstr:tend,1,:,:));
  ubar=squeeze(nc{'ubar'}(tstr:tend,:,:));
  vbar=squeeze(nc{'vbar'}(tstr:tend,:,:));
else
  ui=squeeze(nc{'ubar'}(1,:,:));
  vi=squeeze(nc{'vbar'}(1,:,:));
  mu=squeeze(mean(nc{'ubar'}(tstr:tend,:,:)));
  mv=squeeze(mean(nc{'vbar'}(tstr:tend,:,:)));
  u=nc{'ubar'}(tstr:tend,:,:);
  v=nc{'vbar'}(tstr:tend,:,:);
end


%% 

[~,xn] = min(abs(y(:,1)-60));
xinit = y(xn,1);
if dyn=="3D"
    load /Users/simon/Code/CONFIGS/IB09_S10_newPS/IB09_S10_newPS_avg_Xall_Yall_Tall_Z9_tpas01.mat
    t=out;
    load /Users/simon/Code/CONFIGS/IB09_S10_newPS/IB09_S10_newPS_2_avg_Xall_Yall_Tall_Z9_tpas01.mat
    t=cat(1,t,out);
    load /Users/simon/Code/CONFIGS/IB09_S10_newPS/IB09_S10_newPS_3_avg_Xall_Yall_Tall_Z9_tpas01.mat
    t=cat(1,t,out);
elseif dyn=="2D"
    load /Users/simon/Code/CONFIGS/IB09_DB1tsdiff/IB09_DB1tsdiff_avg_Xall_Yall_Tall_Z9_tpas01.mat
    t=out;
    load /Users/simon/Code/CONFIGS/IB09_DB1tsdiff/IB09_DB1tsdiff_2_avg_Xall_Yall_Tall_Z9_tpas01.mat
    t=cat(1,t,out);
end
t=squeeze(t(:,1,:,:));
[dyr,ix0] = min(abs(x(1,:)+8));
[dyr,ixmax] = min(abs(x(1,:)+356));
time=time-time(1); %spin-up

[dyr,x82] = min(abs(y(:,1)-82-xinit));
[~,tinit] = min(abs(time-2300));
[~,tend] = min(abs(time-19000));
usedtime82=(time(tend)-time(tinit))/(19000-4300);
t82 = squeeze(mean(t(tinit:tend,x82,ixmax:ix0)));
std82 = squeeze(std(t(tinit:tend,x82,ixmax:ix0)));
t82h = t82 + std82./length(tinit:tend);
t82b = t82 - std82./length(tinit:tend);

[dyr,x248] = min(abs(y(:,1)-208-xinit));
[~,tinit] = min(abs(time-1000));
[~,tend] = min(abs(time-23000));
usedtime248=(time(tend)-time(tinit))/(23000-1500);
t248 = squeeze(mean(t(tinit:tend,x248,ixmax:ix0)));
std248 = squeeze(std(t(tinit:tend,x248,ixmax:ix0)));
t248h = t248 + std248/length(tinit:tend);
t248b = t248 - std248/length(tinit:tend);

[dyr,x546] = min(abs(y(:,1)-546-xinit));
[~,tinit] = min(abs(time-5000));
[~,tend] = min(abs(time-24500));
usedtime546=(time(tend)-time(tinit))/(24500-5000);
t546 = squeeze(mean(t(tinit:tend,x546,ixmax:ix0)));
std546 = squeeze(std(t(tinit:tend,x546,ixmax:ix0)));
t546h = t546 + std546/length(tinit:tend);
t546b = t546 - std546/length(tinit:tend);

[dyr,x1069] = min(abs(y(:,1)-1069-xinit));
[~,tinit] = min(abs(time-6000));
[~,tend] = min(abs(time-24500));
usedtime1069=(time(tend)-time(tinit))/(24500-6000);
t2=t(tinit:tend,x1069-0:x1069+0,ixmax:ix0);
t2=squeeze(mean(t2,2));
t1069 = squeeze(mean(t2))';
std1069 = smooth(squeeze(std(t2)));
t1069h = t1069 + std1069/length(tinit:tend);
t1069b = t1069 - std1069/length(tinit:tend);

[dyr,x1662] = min(abs(y(:,1)-1662-xinit));
[~,tinit] = min(abs(time-13000));
[~,tend] = min(abs(time-24500));
t1662 = squeeze(mean(t(tinit:tend,x1662,ixmax:ix0)));
std1662 = squeeze(std(t(tinit:tend,x1662,ixmax:ix0)));
t1662h = t1662 + std1662/length(tinit:tend);
t1662b = t1662 - std1662/length(tinit:tend);

%% Loading data from HR15
load Data/HR16_fig10_D82mod.mat
load Data/HR16_fig10_D82obs.mat
load Data/HR16_fig10_D248mod.mat
load Data/HR16_fig10_D248obs.mat
load Data/HR16_fig10_D1069mod.mat
load Data/HR16_fig10_D1069obs.mat
%% Plotting
figure;
ccroco='blue';cObs='black';cFW='red';
wcroco=2.5;wFW=2;wObs=wFW;

% As for now, only SA1 and SA2 can be plotted (domain length in CROCO=
% 600m)
fig=subplot(3,1,1);
x2 = [x(1,ixmax:ix0), fliplr(x(1,ixmax:ix0))];
inBetween = [t82h', fliplr(t82b')];
fill(x2, inBetween,[.9 .9 .9],'linestyle','none');
hold on
xlim([-400,0]);
plot(x(1,ixmax:ix0),t82'.*min(h(1,ixmax:ix0),2.7)./h(1,ixmax:ix0),'LineWidth',wcroco,'Color',ccroco);
%plot(x(1,ixmax:ix0),DstoD(t82',h(ixmax:ix0)));
plot(HR16_fig10_D82mod(:,1),HR16_fig10_D82mod(:,2),'LineWidth',wFW,'Color',cFW);
plot(H16_fig10_D82obs(:,1),H16_fig10_D82obs(:,2),'LineWidth',wObs,'Color',cObs);
xline(-81,lineStyle="--",lineWidth=2);
ylim([0 100])
legend('','CROCO','funwaveC', 'obs','Location','northwest','Interpreter','latex');
ylabel('$\langle D \rangle$  (ppb)','Interpreter','latex');
title("(a) y = "+string(82)+' m ('+string(round(usedtime82,2)*100)+'\% of time)','Interpreter','latex');

subplot(3,1,2);
inBetween = [t248h', fliplr(t248b')];
fill(x2, inBetween,[.9 .9 .9],'linestyle','none');
hold on;
ylim([0 50])
plot(x(1,ixmax:ix0),t248'.*min(h(1,ixmax:ix0),2.7)./h(1,ixmax:ix0),'LineWidth',wcroco,'Color',ccroco);
plot(HR16_fig10_D248mod(:,1),HR16_fig10_D248mod(:,2),'LineWidth',wFW,'Color',cFW);
plot(HR16_fig10_D248obs(:,1),HR16_fig10_D248obs(:,2),'LineWidth',wObs,'Color',cObs);
xline(-81,lineStyle="--",lineWidth=2);
xlim([-400,0]);
legend('','CROCO','funwaveC', 'obs','Location','northwest','Interpreter','latex');
ylabel('$\langle D \rangle$ (ppb)','Interpreter','latex');
title(strcat("(b) y = ",string(248),' m (',string(round(usedtime248,2)*100),'\% of time)'),'Interpreter','latex');
% 
% % subplot(5,1,3);
% [d,ix] = min(abs(y(:,1)-546-68));
% inBetween = [t546h', fliplr(t546b')];
% fill(x2, inBetween,[.9 .9 .9],'linestyle','none');
% hold on;
% plot(x(1,ixmax:ix0),t546,'Color','red');
% xlim([-350,0]);
% xline(-81,lineStyle="--",lineWidth=2);
% legend('','CROCO mod','Location','northwest','Interpreter','latex');
% ylab = ylabel(' D  (ppb)','Interpreter','latex');
% title(strcat("(c) y = ",string(546),' m (',string(round(usedtime546,2)*100),'\% of time)'),'Interpreter','latex');

subplot(3,1,3);
inBetween = [t1069h', fliplr(t1069b')];
fill(x2, inBetween,[.9 .9 .9],'linestyle','none');
hold on;
plot(x(1,ixmax:ix0),t1069'.*min(h(1,ixmax:ix0),2.7)./h(1,ixmax:ix0),'LineWidth',wcroco,'Color',ccroco);
plot(HR16_fig10_D1069mod(:,1),HR16_fig10_D1069mod(:,2),'LineWidth',wFW,'Color',cFW);
plot(HR16_fig10_D1069obs(:,1),HR16_fig10_D1069obs(:,2),'LineWidth',wObs,'Color',cObs);
xline(-81,lineStyle="--",lineWidth=2);
xlim([-400,0]);
legend('','CROCO','funwaveC', 'obs','Location','northwest','Interpreter','latex');
[d,ix] = min(abs(y(:,1)-1069));
ylab = ylabel(' D  (ppb)','Interpreter','latex');
title(strcat("(c) y = ",string(1069),' m (',string(round(usedtime1069,2)*100),'\% of time)'),'Interpreter','latex');
% 
% subplot(5,1,5);
% [d,ix] = min(abs(y(:,1)-1616));
% inBetween = [t1662h', fliplr(t1662b')];
% fill(x2, inBetween,[.9 .9 .9],'linestyle','none');
% hold on;
% plot(x(1,ixmax:ix0),t1662);
% xlim([-350,0]);
% xline(-81,lineStyle="--",lineWidth=2);
% legend('','mod','Location','northwest','Interpreter','latex');
% ylab = ylabel(' D  (ppb)','Interpreter','latex');
% title(strcat("(d) y = ",string(1662),' m'),'Interpreter','latex');
% xlabel('$x$ (m)','Interpreter','latex');
% saveas(fig,'HR16_fig10.pdf')

print(gcf, '-dpdf', '-bestfit', './profiles_tracer.pdf');

