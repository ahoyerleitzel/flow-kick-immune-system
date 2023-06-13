%% Figures 5 for infectious disease modeling paper
%% load functions

paramset0 = [0.35,20,0.2,0.04,5e-7,2e-5,0.8,7e7,3,1*10^-7,2,1e-4,0.01,12000,0.2,10,3];
fig5 = figure_5_repeatedkickscalculationMay2023;

%% approximate boundary with bisection method
IC = zeros(6,1);
IC(2) = 7e7; 

errtol = 10e-6;
maxtime = 500;
timelist = 2:25;
kickinputs = 15000*ones(size(timelist));

[timelist,kvec] = fig5.curve(timelist,IC,kickinputs,paramset0,maxtime,errtol)
kb = kvec(:,1);
ka = kvec(:,2);
%save('figure_5_boundary_approximation.mat','timelist','ka','kb')
%% alternatively, just load boundary approximation
load('figure_5_boundary_approximation.mat')

%% More detailed boundary calculation
%needs figure_5_boundary_approximation.mat to be loaded
%as well as fig5 = figure_5_repeatedkickscalculationMay2023.m
% Warning: This section is very computationally intensive at the error
% tolerance and grid size used for computation.

IC = zeros(6,1);
IC(2) = 7e7;

errtol = 10e-6; %used 10e-11 for publication
maxtime = 1500;
deltakick = 5; %grid step in kick %used deltak = 1 for publication
deltatau = 1; %grid step in tau %used deltatau = 0.1 for publication

taulist = [2,5:5:30];
klist = [ka(1),ka(4:5:end)',0];
for j1 = 5:length(taulist)-1
filename = ['figure5_tau',num2str(taulist(j1)),'to',num2str(taulist(j1+1)),'.mat']
tausteps = taulist(j1):deltatau:taulist(j1+1); %delta tau steps
slope = (klist(j1+1)-klist(j1))/(taulist(j1+1)-taulist(j1)); %slope of linear approximiation to boundary
parfor k1 = 1:length(tausteps)
    tic
    tauvalue = tausteps(k1)
    counter = 1;
    thresh = 1.01;
    excursionkick = [];
    postkickstates = [];
    % set how far above or below boundary approximation
    bottomkick = max((tauvalue-taulist(j1))*slope + klist(j1) - 1000,10);
    topkick = min((tauvalue-taulist(j1))*slope + klist(j1) + 5000,16500);
    for kick = bottomkick:deltakick:topkick
        [exc,ICnew] = fig5.kickflowexcursion(IC, tauvalue, kick, thresh, paramset0, ceil(maxtime/tauvalue),errtol);
        excursionkick(counter,1) = kick;
        excursionkick(counter,2) = exc;
        counter = counter + 1;
    end
    val = excursionkick(excursionkick(:,2)==1,1);
    excursioninfo(k1).excursionvalues = val;
    val2 = excursionkick(excursionkick(:,2)==0,1);
    excursioninfo(k1).noexcursionvalues = val2;
    excursioninfo(k1).excursionkick = excursionkick;
    excursioninfo(k1).tau = tauvalue;
    toc
end
save(filename,'excursioninfo')
end
%% make figure
% load .mat files used for publication figure
files = dir('figure5_*refined.mat')
color1 = [0, 0.4470, 0.7410];
min_k = [];
taulist = [];
scale = 1e5;
for i = 1:length(files)
    load(files(i).name)
    for k1 = 1:length(excursioninfo)
        val1 = excursioninfo(k1).excursionvalues;
        val2 = excursioninfo(k1).noexcursionvalues;
        tau = excursioninfo(k1).tau;
        hold on
        index = find(val2 > 100);
        % not plotting the bottom for a nicer image.
        plot(tau+0*val2(index),val2(index)/scale,'.', 'Color', color1,'MarkerSize',4)
        taulist(end+1) = tau;
        min_k(end+1) = min(val2);
    end
end
% fill in area below boundary
area([1.7 1.7 taulist 30 1.97], [10 .095 (min_k+10) 10 10]/scale,'FaceColor',color1,'EdgeColor','none')

%% add ticks and box, save as eps
axis([0 30 0 0.18])
xticks([0:5:30])
yticks([0:0.03:0.18])
ax.XAxis.FontSize = 6;
ax.YAxis.FontSize = 6;
grid on
xlabel("flow time \tau", 'FontSize', 8)
ylabel("kick size x1e5",'FontSize',8)
set(gca,'LooseInset',get(gca,'TightInset'));
hold on
rectangle('Position',[7 12000/scale 1 1000/scale],'FaceColor',"#EDB120")


eps_save(1,'fig5.eps')
%%
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