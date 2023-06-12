%% Figure 3 - single excursion
% Use an initial virus level of 4000
% Kick once after time tau
% Look for an excursion (or not)

function funcs = figure_3g_onekickexcursioncalculationMay2023
% Call this function first
% This gives access to all the functions in this file from the command line.
% Call with syntax of funcs = Figure3_onekickexcursioncalculationMay2023;
funcs.flow = @flow;
funcs.checkforexcursion = @checkexcursion;
funcs.kickflowexcursion = @kickflowexcursion;
funcs.findET = @findexcursionthreshold;
funcs.plot = @plotflowkick;
funcs.curve = @thresholdcurve;
end

% Original parameter set
% paramset0 = [0.35,20,0.2,0.04,5e-7,2e-5,0.8,7e7,3,1*10^-7,2,1e-4,0.01,12000,0.2,10,3];

function [tvec,kvec] = thresholdcurve(timevector,IC,kickinputs,numkicks,paramset)
% This function finds a threshold value for vector of flow times (timevector)
% The threshold value is kick size that does not cause excursion with a 
%   nearby kick size that does cause it
% 
tvec=[];
kvec=[];
for k1=1:length(timevector)
    tau = timevector(k1);
    kick0 = kickinputs(k1);
    [belowthresh,abovethresh] = findexcursionthreshold(tau,IC,kick0,numkicks,paramset);
    tvec = [tvec,tau]
    kvec = [kvec;[belowthresh,abovethresh]];
end
end

function plotflowkick(IC,kickv,tau,params)
% Plot flows and kicks.  Use plots to check the calculations.
allt = [0];
allY = zeros(1,6);

% Initial run:
[t1,Y1] = flow(tau, IC, params);
allt(end+1:end+length(t1)) = t1 + allt(end);
allY(end+1:end+length(Y1),:) = Y1;

% Kick once
IC = Y1(end,:);
IC(1) = IC(1) + kickv;

% Leave time for the excursion to develop after the kick
flowtime = max(tau,40);
[t1,Y1] = flow(flowtime, IC, params);
allt(end+1:end+length(t1)) = t1 + allt(end);
allY(end+1:end+length(Y1),:) = Y1;

% Plot
plot(allt, allY(:,1))
end


function [below,above] = findexcursionthreshold(tau,IC,kickinit,numkicks,paramset)
% Initial condition of the first flow is fixed for these calculations
% The kick size is varied to find the excursion.
kick = kickinit;
thresh = 0.5;
exclist = kickflowexcursion(IC,tau,kick,thresh,paramset);
% Loop until we do see an excursion
counter = 1;
while sum(exclist)==0 && counter < 20
    kick = kick*2;
    exclist = kickflowexcursion(IC,tau,kick,thresh,paramset);
    counter = counter + 1;
end
if sum(exclist)==0
    disp('no excursion')
    aboveexcursion = 0;
    belowexcursion = 0;
else
    aboveexcursion = kick;
    belowexcursion = 0;
end
counter = 1;
while (aboveexcursion-belowexcursion)>1 && counter < 15
    kick = (aboveexcursion + belowexcursion)/2;
    exclist = kickflowexcursion(IC,tau,kick,thresh,paramset);
    if sum(exclist)>0
        aboveexcursion = kick;
    else
        belowexcursion = kick;
    end
    counter = counter + 1;
end
below = belowexcursion;
above = aboveexcursion;
end


function exc = checkexcursion(t,Y1,kick,thresholdfactor)
% Takes in a timeseries of system state values (t,Y1), the initial kick, 
%   and a threshold factor.
% The virus values are stored in Y1(:,1)
% Check for virus ever being above the "excursion threshold"

virus = Y1(:,1);
[~,index] = max(diff(virus)./diff(t));
% Check that after the positive difference the virus reaches the threshold:
excursion = max(virus(index:end,1));
% ensure the threshold is > 0 even if the kick is 0.
thresh = max(kick*thresholdfactor,1);
if excursion > thresh
        exc = 1;
    else
        exc = 0;
    end
end


function exc1 = kickflowexcursion(IC, tau, kick, thresholdfactor,paramset)
% Initial run:
[~,Y1] = flow(tau, IC, paramset);
% Kick once
IC = Y1(end,:);
IC(1) = IC(1) + kick;
% Leave time for the excursion to develop after the kick
flowtime = max(tau,40);
[t,Y1] = flow(flowtime, IC, paramset);
exc1 = checkexcursion(t,Y1,IC(1),thresholdfactor);
end


function dydtval = dydt(t,Y,params)
% Parameters are provided in params.
% Order is
% [p,c,mu,muprime,beta,betaprime,g,Ct,delta,q,d,m1,m2,m3,r,Vm,kappa]
p = params(1);
c = params(2);
mu =params(3);
muprime = params(4);
beta =params(5);
betaprime = params(6);
g = params(7);
Ct = params(8);
delta = params(9);
q = params(10);
d = params(11);
m1 = params(12);
m2 = params(13);
m3 = params(14);
r = params(15);
Vm = params(16);
kappa = params(17);

% variable order is 
%  [virus,target cells,infected cells,interferon,B cells,antibodies]
dydtval = [p*Y(3)-c*Y(1)-mu*Y(1)*Y(6)-beta*Y(1)*Y(2)*Y(1)/(Vm+Y(1));...
    g*(Y(2))*(1-(Y(2)+Y(3))/(Ct))-betaprime*Y(1)*Y(2)*Y(1)/(Vm+Y(1));...
    betaprime*Y(1)*Y(2)*Y(1)/(Vm+Y(1))-delta*Y(3)-kappa*Y(3)*Y(4);...
    q*Y(3)-d*Y(4);...
    m1*Y(1)*(1-Y(5))-m2*Y(5);...
    m3*Y(5)-r*Y(6)-muprime*Y(1)*Y(6)];
end

function [ts,YS] = flow(tau,IC,params)
options = odeset('AbsTol',1e-6,'RelTol',1e-6);
[ts,YS] = ode15s(@(t,y) dydt(t,y,params),[0,tau],IC,options);
end