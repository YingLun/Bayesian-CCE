% y(t) = A * x(t) + e(t)       e(t) ~ N(0,R)
% x(t) = C x(t-1) + u(t)       u(t) ~ N(0,Q)
%% Housekeeping
clear('all'); close('all'); clc;
main_path = 'C:\Users\Pooyan\Dropbox\Teaching\WS2014\Bayesiann Econometrics\MatlabExcercise\KalmanFilterSimpleExample';
cd(main_path)
%% Simulate Data
T       = 200;
N       = 1;
K       = 1;
yL      = zeros(N,T);
yH      = zeros(N,T);
x       = zeros(K,T);
% sigma_e = [.1 0;0 .08];
sigma_eL = 1;
sigma_eH = 5;
sigma_u  = 1;
eL       = mvnrnd(0,sigma_eL,T);
eH       = mvnrnd(0,sigma_eH,T);
u        = mvnrnd(0,sigma_u,T);
A        = 1;%[.85 -.75];
B        = .9;

for t = 2:T
    x(:,t)  = B * x(:,t-1) + u(t);
    yL(:,t) = A * x(:,t) + eL(t);
    yH(:,t) = A * x(:,t) + eH(t);
end

%% Kalman Filter
% =========================================================================
[kf_outputL] = kf_class(yL,A,B,sigma_eL,sigma_u) ;
[kf_outputH] = kf_class(yH,A,B,sigma_eH,sigma_u) ;
% kf_output.Kt_gain  = Kt_gain;     Kalman Gain
% kf_output.Xtt1     = Xtt1;        Predicted State
% kf_output.Xtt      = Xtt;         Updated State
% kf_output.P_tt     = P_tt;        Cov(Predicted State)
% kf_output.P_tt1    = P_tt1;       Cov(Updated State)


%% Kalman Smoother
% =========================================================================
[ks_outputL] = ks_class(kf_outputL,B);
[ks_outputH] = ks_class(kf_outputH,B);
% ks_output.XtT      = Xtt;         Smoothed State
% ks_output.P_tT     = P_tt;        Cov(Smoothed State)


min_val = floor(10*min([yL(:);yH(:)]))/10;
max_val = ceil(10*max([yL(:);yH(:)]))/10;

min_val1 = floor(10*min(x))/10;
max_val1 = ceil(10*max(x))/10;


%% Figure 1: Low Case
Kt_gain = kf_outputL.Kt_gain;
Xtt1    = kf_outputL.Xtt1;
Xtt     = kf_outputL.Xtt;
XtT     = ks_outputL.XtT;

FS = 12; LW = 3;

figure('Name','DataStateLow');orient('Landscape')
subplot(2,1,1);hold('on')
plot(1:T,yL,'r-','LineWidth',2);
plot(1:T,x,'color',[.8 .8 .8]-.3,'Linewidth',4)
grid;% axis('tight');
title('Unobserved State (Low)','Fontsize',FS)
legend('Data: Y','State: X_{t}','Location','NorthWest')
grid('on');set(gca,'FontSize',FS,'FontWeight','Bold','YLim',[min_val max_val])

subplot(2,1,2);hold('on')
plot(x(1,1:end),'color',[.8 .8 .8]-.3,'Linewidth',4)
plot(Xtt(1,1:end),'-b','Linewidth',LW)
legend('State: X_{t}','Updated: X_{t|t}','Location','NorthWest')
grid('on');set(gca,'FontSize',FS,'FontWeight','Bold','YLim',[min_val max_val])
saveas(gcf,'DataStateLow','pdf')

figure('Name','StateLow');orient('Landscape')
hold on
plot(x(1,1:end),'color',[.8 .8 .8]-.3,'Linewidth',3)
plot(Xtt1(1,1:end),'-*r','Linewidth',LW)
plot(Xtt(1,1:end),'-b','Linewidth',LW)
plot(XtT(1,1:end),'--k','Linewidth',LW+1)
title('Unobserved State (Low)','Fontsize',FS)
legend('True: X','Predicted: X_{t|t-1}','Updated: X_{t|t}','Smoothed; X_{t|T}','Location','NorthEast')
grid('on')
set(gca,'FontSize',FS,'FontWeight','Bold','LineWidth',LW,'YLim',[min_val1 max_val1])
saveas(gcf,'StateLow','pdf')

%% Low case Figure 1
figure('Name','StateLow1');orient('Landscape')
hold on
plot(x(1,1:end),'color',[.8 0 0],'Linewidth',3)
title('Unobserved State (Low)','Fontsize',FS)
legend('True: X','Location','NorthEast')
grid('on')
set(gca,'FontSize',FS,'YLim',[min_val1 max_val1])
saveas(gcf,'StateLow1','pdf')
close('all')

figure('Name','StateLow2');orient('Landscape')
hold on
plot(x(1,1:end),'color',[.8 0 0],'Linewidth',3)
plot(Xtt1(1,1:end),'-','Color',[.8 .8 .8]-.2,'Linewidth',LW)
title('Unobserved State (Low)','Fontsize',FS)
legend('True: X','Predicted: X_{t|t-1}','Location','NorthEast')
grid('on')
set(gca,'FontSize',FS,'YLim',[min_val1 max_val1])
saveas(gcf,'StateLow2','pdf')
close('all')

figure('Name','StateLow3');orient('Landscape')
hold on
plot(x(1,1:end),'color',[.8 0 0],'Linewidth',3)
plot(Xtt1(1,1:end),'-','Color',[.8 .8 .8]-.2,'Linewidth',LW)
plot(Xtt(1,1:end),'-b','Color',[.8 .8 .8]-.4,'Linewidth',LW)
title('Unobserved State (Low)','Fontsize',FS)
legend('True: X','Predicted: X_{t|t-1}','Updated: X_{t|t}','Location','NorthEast')
grid('on')
set(gca,'FontSize',FS,'YLim',[min_val1 max_val1])
saveas(gcf,'StateLow3','pdf')
close('all')

figure('Name','StateLow4');orient('Landscape')
hold on
plot(x(1,1:end),'color',[.8 0 0],'Linewidth',3)
plot(Xtt1(1,1:end),'-','Color',[.8 .8 .8]-.2,'Linewidth',LW)
plot(Xtt(1,1:end),'-','Color',[.8 .8 .8]-.4,'Linewidth',LW)
plot(XtT(1,1:end),'-','Color',[.8 .8 .8]-.6,'Linewidth',LW+1)
title('Unobserved State (Low)','Fontsize',FS)
legend('True: X','Predicted: X_{t|t-1}','Updated: X_{t|t}','Smoothed; X_{t|T}','Location','NorthEast')
grid('on')
set(gca,'FontSize',FS,'YLim',[min_val1 max_val1])
saveas(gcf,'StateLow4','pdf')
close('all')

figure('Name','StateLow5');orient('Landscape')
hold on
plot(x(1,1:end),'color',[.8 0 0],'Linewidth',3)
plot(Xtt1(1,1:end),'-','Color',[.8 .8 .8]-.2,'Linewidth',LW)
plot(Xtt(1,1:end),'-','Color',[.8 .8 .8]-.4,'Linewidth',LW)
plot(XtT(1,1:end),'-','Color',[.8 .8 .8]-.6,'Linewidth',LW+1)
title('Unobserved State (Low)','Fontsize',FS)
legend('True: X','Predicted: X_{t|t-1}','Updated: X_{t|t}','Smoothed; X_{t|T}','Location','NorthEast')
grid('on')
set(gca,'FontSize',FS,'XLim',[1 50],'YLim',[min_val1 max_val1])
saveas(gcf,'StateLow5','pdf')
close('all')

%% ========================================================================
%% Figure 1: Low High
Kt_gain = kf_outputH.Kt_gain;
Xtt1    = kf_outputH.Xtt1;
Xtt     = kf_outputH.Xtt;
XtT     = ks_outputH.XtT;

FS = 12; LW = 3;

figure('Name','DataStateHigh');orient('Landscape')
subplot(2,1,1);hold('on')
plot(1:T,yH,'r-','LineWidth',2);
plot(1:T,x,'color',[.8 .8 .8]-.3,'Linewidth',4)
grid;% axis('tight');
title('Unobserved State (High)','Fontsize',FS)
legend('Data: Y','State: X_{t}','Location','NorthWest')
grid('on');set(gca,'FontSize',FS,'FontWeight','Bold','YLim',[min_val max_val])
subplot(2,1,2);hold('on')
plot(x(1,1:end),'color',[.8 .8 .8]-.3,'Linewidth',4)
plot(Xtt(1,1:end),'-b','Linewidth',LW)
legend('State: X_{t}','Updated: X_{t|t}','Location','NorthWest')
grid('on');set(gca,'FontSize',FS,'FontWeight','Bold','YLim',[min_val max_val])
saveas(gcf,'DataStateHigh','pdf')

figure('Name','StateHigh');orient('Landscape')
hold on
plot(x(1,1:end),'color',[.8 .8 .8]-.3,'Linewidth',3)
plot(Xtt1(1,1:end),'-*r','Linewidth',LW)
plot(Xtt(1,1:end),'-b','Linewidth',LW)
plot(XtT(1,1:end),'--k','Linewidth',LW+1)
title('Unobserved State (High)','Fontsize',FS)
legend('True: X','Predicted: X_{t|t-1}','Updated: X_{t|t}','Smoothed; X_{t|T}','Location','NorthEast')
grid('on')
set(gca,'FontSize',FS,'FontWeight','Bold','LineWidth',LW,'YLim',[min_val1 max_val1])
saveas(gcf,'StateHigh','pdf')


%% Low case Figure 1
figure('Name','StateHigh1');orient('Landscape')
hold on
plot(x(1,1:end),'color',[0 0 .7],'Linewidth',3)
title('Unobserved State (High)','Fontsize',FS)
legend('True: X','Location','NorthEast')
grid('on')
set(gca,'FontSize',FS,'YLim',[min_val1 max_val1])
saveas(gcf,'StateHigh1','pdf')
close('all')

figure('Name','StateHigh2');orient('Landscape')
hold on
plot(x(1,1:end),'color',[0 0 .7],'Linewidth',3)
plot(Xtt1(1,1:end),'-','Color',[.8 .8 .8]-.2,'Linewidth',LW)
title('Unobserved State (High)','Fontsize',FS)
legend('True: X','Predicted: X_{t|t-1}','Location','NorthEast')
grid('on')
set(gca,'FontSize',FS,'YLim',[min_val1 max_val1])
saveas(gcf,'StateHigh2','pdf')
close('all')

figure('Name','StateHigh3');orient('Landscape')
hold on
plot(x(1,1:end),'color',[0 0 .7],'Linewidth',3)
plot(Xtt1(1,1:end),'-','Color',[.8 .8 .8]-.2,'Linewidth',LW)
plot(Xtt(1,1:end),'-b','Color',[.8 .8 .8]-.4,'Linewidth',LW)
title('Unobserved State (High)','Fontsize',FS)
legend('True: X','Predicted: X_{t|t-1}','Updated: X_{t|t}','Location','NorthEast')
grid('on')
set(gca,'FontSize',FS,'YLim',[min_val1 max_val1])
saveas(gcf,'StateHigh3','pdf')
close('all')

figure('Name','StateHigh4');orient('Landscape')
hold on
plot(x(1,1:end),'color',[0 0 .7],'Linewidth',3)
plot(Xtt1(1,1:end),'-','Color',[.8 .8 .8]-.2,'Linewidth',LW)
plot(Xtt(1,1:end),'-','Color',[.8 .8 .8]-.4,'Linewidth',LW)
plot(XtT(1,1:end),'-','Color',[.8 .8 .8]-.6,'Linewidth',LW+1)
title('Unobserved State (High)','Fontsize',FS)
legend('True: X','Predicted: X_{t|t-1}','Updated: X_{t|t}','Smoothed; X_{t|T}','Location','NorthEast')
grid('on')
set(gca,'FontSize',FS,'YLim',[min_val1 max_val1])
saveas(gcf,'StateHigh4','pdf')
close('all')

figure('Name','StateHigh5');orient('Landscape')
hold on
plot(x(1,1:end),'color',[0 0 .7],'Linewidth',3)
plot(Xtt1(1,1:end),'-','Color',[.8 .8 .8]-.2,'Linewidth',LW)
plot(Xtt(1,1:end),'-','Color',[.8 .8 .8]-.4,'Linewidth',LW)
plot(XtT(1,1:end),'-','Color',[.8 .8 .8]-.6,'Linewidth',LW+1)
title('Unobserved State (High)','Fontsize',FS)
legend('True: X','Predicted: X_{t|t-1}','Updated: X_{t|t}','Smoothed; X_{t|T}','Location','NorthEast')
grid('on')
set(gca,'FontSize',FS,'XLim',[1 50],'YLim',[min_val1 max_val1])
saveas(gcf,'StateHigh5','pdf')
close('all')


