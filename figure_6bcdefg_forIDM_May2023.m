%% Figure 6bcdefg
% note this code uses the "old" version of the underlying ODEs with seven
% variables. The resistance cells (the 4th variable) is set to 0 and left
% out of the model.
Ct = 7e7;
scale = [10^5, Ct, Ct, 1, 1, 1, 5*10^4];
IC0 = zeros(1,7);
IC0(2) = Ct;
titles = ["Viral Load x1e5", "Target Cells x7e7", "Infected Cells x7e7", "Interferon", "B-cells", "Antibodies x5e4"];

%% uniform box (7,8)x(12000,13000)
maxsims = 500;
mintau=7;
excursionlist = zeros(1,500);
for j = 1:maxsims
IC0(1) = 15000;
IC = IC0;
ts1 = [];
YS1 = [];
taulist = [];
kicklist = [];
for i = 1:ceil(600/mintau)
   tau = 7+1*rand;
   taulist = [taulist, tau];
   [t,Y] = flow(tau,IC);
   if i==1
       ts1(end+1:end+length(t)) = t;
   else
       ts1(end+1:end+length(t)) = t+ts1(end);
   end
   YS1(end+1:end+length(t),:) = Y;
   kicksize = 12000 + 1000*rand();
   kicklist = [kicklist, kicksize];
   ykick = kick(kicksize, Y(end,:));
   IC = ykick;
   exc = checkexcursion(ts1,YS1);
   excursionlist(j) = exc;
end
sim(j).exc = exc;
sim(j).ts = ts1;
sim(j).YS = YS1;
sim(j).taulist = taulist;
sim(j).kicklist = kicklist;
end

excursions = [];
for j=1:maxsims
    excursions = [excursions, sim(j).exc];
end
sum(excursions)


%% plot
figure(1) 
tlength = 600;
tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
label = ['B','C','D','E','F','G'];
for i=1:3
    nexttile
    for j=1:50
        plot(sim(j).ts,sim(j).YS(:,i)/scale(i),'LineWidth',1)
        hold on
    end
    if i == 1
        axis([0 tlength 0 3])
        yticks(0:1.5:3)
        yticklabels({"0","","3"})
        text(0,3.5,label(i),'EdgeColor','k','fontsize',8,'Margin',0.5)
    else
        axis([ 0 tlength 0 1.1])
        yticks(0:0.5:1)
        yticklabels({"0","","1"})
        text(0,1.3,label(i),'EdgeColor','k','fontsize',8,'Margin',0.5)
    end
    xticks(0:100:600)
    set(gca,'FontSize',6)
    title(titles(i),'fontsize',8)
    grid on
end

for i=4:6
    nexttile
    for j=1:50
        plot(sim(j).ts,sim(j).YS(:,i+1)/scale(i+1),'LineWidth',1)
        hold on
    end
    axis([ 0 tlength 0 1.1])
    xticks(0:100:600)
    yticks(0:0.5:1)
    yticklabels({"0","","1"})
    set(gca,'FontSize',6)
    title(titles(i),'fontsize',8)
    text(0,1.3,label(i),'EdgeColor','k','fontsize',8,'Margin',0.5)
    grid on
end

eps_save(1,'fig6bcdefg.eps')


%%
function exc = checkexcursion(ts,YS)
% takes in a timeseries of system state values, the initial kick, and a
% threshold factor
% virus is in the first column of the system values
% check for virus ever being above the "excursion threshold"
    tequals10 = find(ts >= 10,1);
    excursion = max(YS(tequals10:end,1)); % this is the max virus after t=10.
    if excursion > 2*10^5
        exc = 1;
    else
        exc = 0;
    end
end


%%
function [YS] = kick(kicksize,IC)
    YS = IC;
    YS(1) = YS(1)+kicksize;
end

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
phi = 0; %0.14 or 0
rho = 0; %0.05 or 0
xi = 0; %0.1 or 0
s = 0; %1 or 0
kappa = 3; %3 or 0
q = 1*10^-7; %q0*10e-7
d = 2;
m1 = 1e-4;
m2 = 0.01;
m3 = 12000;
r = 0.2;
Vm = 10;

%[v,T,i,r,f,B,a]
f = @(t,Y) [p*Y(3)/(1+s*Y(5))-c*Y(1)-mu*Y(1)*Y(7)-beta*Y(1)*Y(2)*Y(1)/(Vm+Y(1));...
    g*(Y(2)+Y(4))*(1-(Y(2)+Y(4)+Y(3))/(Ct))-betaprime*Y(1)*Y(2)*Y(1)/(Vm+Y(1))-rho*Y(4)-phi*Y(5)*Y(2);...
    betaprime*Y(1)*Y(2)*Y(1)/(Vm+Y(1))-delta*Y(3)-kappa*Y(3)*Y(5);...
    phi*Y(5)*Y(2)-rho*Y(4)-xi*Y(4);...
    q*Y(3)-d*Y(5);...
    m1*Y(1)*(1-Y(6))-m2*Y(6);...
    m3*Y(6)-r*Y(7)-muprime*Y(1)*Y(7)];

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
