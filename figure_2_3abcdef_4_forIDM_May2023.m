%% Figures 2, 3ABCDEF, 4 for infectious disease modeling paper

%% Initial conditions and plot scaling factors.
Ct = 7e7;
scale = [10^5, Ct, Ct, 1, 1, 5*10^4];
IC0 = zeros(1,6);
IC0(2) = Ct;
IC0(1) = 1;
titles = ["Viral Load 1xe5", "Target Cells x7e7", "Infected Cells x7e7", "Interferon", "B-cells", "Antibodies x5e4"];
label1 = ['A','B','C','D','E','F'];
label2 = ['G', 'H', 'I', 'J', 'K', 'L'];

%% Figure 2 - Set virus to 15000 and integrate forward in time
initialvirus = 15000;
IC0(1) = initialvirus;

integrationtime = 650;
[ts,YS] = flow(integrationtime,IC0);

% Figure 2 - Plot six time series

figure(1)
tlength = 10;
tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
for k1=1:6
    nexttile
    plot(ts,YS(:,k1)/scale(k1),'LineWidth',1)
if k1 == 1
    axis([0 tlength 0 3])
    yticks(0:1.5:3)
    yticklabels({"0","","3"})
    text(0,3.5,label1(k1),'EdgeColor','k','fontsize',8,'Margin',0.5)
else
    axis([ 0 tlength 0 1.1])
    yticks(0:0.5:1)
    yticklabels({"0","","1"})
    text(0,1.3,label1(k1),'EdgeColor','k','fontsize',8,'Margin',0.5)
end
    xticks(0:2:10)
    set(gca,'FontSize',6)
    title(titles(k1),'fontsize',8)
    
    grid on
end
eps_save(1,'fig2-flow.eps')

%% Figure 3 - Add a second kick that gives an excursion

tau = 283;
kicksize = 15000;
IC0(1) = kicksize;
IC = IC0;
ts1 = [];
YS1 = [];
for k1 = 1:2
   [t,Y] = flow(tau,IC);
   ts1(end+1:end+length(t)) = t + tau*(k1-1);
   YS1(end+1:end+length(t),:) = Y;
   ICnew = Y(end,:);
   ICnew(1) = ICnew(1) + kicksize;
   IC = ICnew;
end
ts_high = ts1;
YS_high = YS1;

% slightly smaller kick at same time that does not give an excursion
ts = ts1;
YS = YS1;
tau = 283;
kicksize = 14800;
IC0(1) = kicksize;
IC = IC0;
ts1 = [];
YS1 = [];
for k1 = 1:2
   [t,Y] = flow(tau,IC);
   ts1(end+1:end+length(t)) = t + tau*(k1-1);
   YS1(end+1:end+length(t),:) = Y;
   ICnew = Y(end,:);
   ICnew(1) = ICnew(1) + kicksize;
   IC = ICnew;
end
ts_low = ts1;
YS_low = YS1;

%%
% Figure 3 - Part ABCDEF.  Plot kick-flow figures above and below excursion
tlength = 400;
label3 = ['A','D','B','E','C','F']
figure(1)
tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');

% Plot time-series
for k1=[1,5,6]
    nexttile
    plot(ts_low,YS_low(:,k1)/scale(k1),'LineWidth',1)
    title(titles(k1),'fontsize',8)
    nexttile
    plot(ts_high,YS_high(:,k1)/scale(k1),'LineWidth',1)
    title(titles(k1),'fontsize',8)
end


% Format and label
for k1 = 1:6
    nexttile(k1)
    if k1 == 1 || k1 == 2
        axis([0 tlength 0 3])
        text(0,3.5,label3(k1),'EdgeColor','k','fontsize',8,'Margin',0.5)
        yticks(0:1.5:3)
        yticklabels({"0","","3"})
    else
        axis([ 0 tlength 0 1.1])
        text(0,1.3,label3(k1),'EdgeColor','k','fontsize',8,'Margin',0.5)
        yticks(0:0.5:1)
        yticklabels({"0","","1"})
    end
    set(gca,'FontSize',6)
    xticks([0,100,200,283,300,400])
    xticklabels({"0","100","200","","300","400"})
    grid on
end

% Mark second exposure
for k1 = 1:2
    nexttile(k1)
    hold on
    hxl = xline(283,'r',{'Second', 'Exposure'},'LabelHorizontalAlignment','left');
    hxl.FontSize = 4.5;
end
for k1 = 3:6
    nexttile(k1)
    hold on
    hxl = xline(283,'r');
end

% Save
eps_save(1,'fig3abcdef.eps')

%% Figure 4ABCDEF - no excursion.  Computational runs
% Run to time 400 with kicks every 7 days
% Use a kick size that is below the excursion threshold.

tau = 7;
kicksize = 15000;
IC0(1) = kicksize;
IC = IC0;
ts1 = [];
YS1 = [];
for k1 = 1:ceil(400/tau)
   [t,Y] = flow(tau,IC);
   ts1(end+1:end+length(t)) = t+tau*(k1-1);
   YS1(end+1:end+length(t),:) = Y;
   ICnew = Y(end,:);
   ICnew(1) = ICnew(1) + kicksize;
   IC = ICnew;
end

tsNE = ts1;
YSNE = YS1;

%% Figure 4ABCDEF - plot
figure(1)
tlength = 400;
tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
for k1=1:6
    nexttile
    plot(tsNE,YSNE(:,k1)/scale(k1),'LineWidth',1)
    if k1 ==1
        axis([ 0 tlength 0 3])
        yticks(0:1.5:3)
        yticklabels({"0","","3"})
        text(0,3.5,label1(k1),'EdgeColor','k','fontsize',8,'Margin',0.5)
    else
        axis([ 0 tlength 0 1.1])
        yticks(0:0.5:1)
        yticklabels({"0","","1"})
        text(0,1.3,label1(k1),'EdgeColor','k','fontsize',8,'Margin',0.5)
    end
    xticks(0:100:400)
    set(gca,'FontSize',6)
    title(titles(k1),'fontsize',8)
    grid on
end

eps_save(1,'fig4A-protection.eps')


%% Figure 4GHIJKL - excursion.  Computational runs
% Run to time 400 with kicks every 7 days
% Use a kick size that is above the excursion threshold.

tau = 7;
kicksize = 16244;
IC0(1) = kicksize;
IC = IC0;
ts1 = [];
YS1 = [];
for k1 = 1:ceil(400/tau)
   [t,Y] = flow(tau,IC);
   ts1(end+1:end+length(t)) = t+tau*(k1-1);
   YS1(end+1:end+length(t),:) = Y;
   ICnew = Y(end,:);
   ICnew(1) = ICnew(1) + kicksize;
   IC = ICnew;
end

tsE = ts1;
YSE = YS1;

%% Figure 4GHIJKL - excursion.  Plot.
figure(1)
tlength = 400;
tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
for k1=1:6
    nexttile
    plot(ts1,YS1(:,k1)/scale(k1),'LineWidth',1)
    if k1 ==1
        axis([ 0 tlength 0 3])
        yticks(0:1.5:3)
        yticklabels({"0","","3"})
        text(0,3.5,label1(k1),'EdgeColor','k','fontsize',8,'Margin',0.5)
    else
        axis([ 0 tlength 0 1.1])
        yticks(0:0.5:1)
        yticklabels({"0","","1"})
        text(0,1.3,label1(k1),'EdgeColor','k','fontsize',8,'Margin',0.5)
    end
    xticks(0:100:400)
    set(gca,'FontSize',6)
    title(titles(k1),'fontsize',8)
    grid on
end

eps_save(1,'fig4B-excursion.eps')



%%
function [ts,YS] = flow(tau,IC)
%This are the parameters for Model 3.
p = 0.35;
c = 20;
mu =0.2;
muprime = 0.04;
beta =5e-7;
betaprime = 2e-5;
g = 0.8;
Ct = 7e7;
delta = 3;
kappa = 3; %3 or 0
q = 1*10^-7; %q0*10e-7
d = 2;
m1 = 1e-4;
m2 = 0.01;
m3 = 12000;
r = 0.2;
Vm = 10;

%[v,T,i,f,B,a]
f = @(t,Y) [p*Y(3)-c*Y(1)-mu*Y(1)*Y(6)-beta*Y(1)*Y(2)*Y(1)/(Vm+Y(1));...
    g*(Y(2))*(1-(Y(2)+Y(3))/(Ct))-betaprime*Y(1)*Y(2)*Y(1)/(Vm+Y(1));...
    betaprime*Y(1)*Y(2)*Y(1)/(Vm+Y(1))-delta*Y(3)-kappa*Y(3)*Y(4);...
    q*Y(3)-d*Y(4);...
    m1*Y(1)*(1-Y(5))-m2*Y(5);...
    m3*Y(5)-r*Y(6)-muprime*Y(1)*Y(6)];

options = odeset('AbsTol',1e-6,'RelTol',1e-6);
[ts,YS] = ode15s(f,[0,tau],IC,options);
end

function y = eps_save(fig_number,filename)
figure(fig_number)

set(gcf,'PaperUnits','inches');
oldsizes = get(gcf,'PaperPosition');
% This returns [x y width height]
newwidth = 3.2;
newheight = oldsizes(4)/oldsizes(3)*newwidth;
set(gcf,'PaperPosition',[0 0 newwidth newheight]);
print('-opengl',filename,'-depsc','-r300')
end