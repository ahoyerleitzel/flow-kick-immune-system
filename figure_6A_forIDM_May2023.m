%% Figure 6 stochastic kicks 
% note this code uses the "old" version of the underlying ODEs with seven
% variables. The resistance cells (the 4th variable) is set to 0 and left
% out of the model.

Ct = 7e7;
scale = [1*10^5, Ct, Ct, 1, 1, 1, 5*10^4];
IC0 = zeros(1,7);
IC0(2) = Ct;
titles = ["Viral Load x1e5", "Target Cells x7e7", "Infected Cells x7e7", "Interferon", "B-cells", "Antibodies x5e4"];


%% 6A - one simulation just viral load
%load previously done simulation
%over all basic structure of plot
load('figure_6_stochastic_example.mat')
ax=axes();
plot(ts,YS(:,1),'LineWidth',2)
axis([100 150 11500 13500])
grid on
xticks(100:5:150)
xticklabels({"100","","","","","","","","","","150"})
yticks(11500:500:13500)
yticklabels({"11500","12000","12500","13000","13500"})
xlabel('time','Fontsize',6)
ylabel('kicksize','Fontsize',6)
set(gca,'FontSize',6)
text(100.5,13250,'A','EdgeColor','k','fontsize',8,'Margin',0.5)
hold on
%% add red ticks at kick times
myTick = [107.07,114.09,121.13,128.61,136.14,143.63];

%ax2
ax2 = axes();
axis([100 150 11500 13500])
yticks(11500:500:13500)
yticklabels({"11500","12000","12500","13000","13500"})

ylabel('kicksize','Fontsize',6)
set(gca,'FontSize',6)
set(ax2, 'XTick', myTick);
set(ax2, 'XColor', 'r');
xlabel('time','Fontsize',6,'Color','k')
set(ax2, 'FontSize',6);


%% cover up red axis with a black axis
%ax3
ax3 = axes();
axis([100 150 11500 13500])
xticks(100:5:150)
xticklabels({"100","","","","","","","","","","150"})
yticks(11500:500:13500)
yticklabels({"11500","12000","12500","13000","13500"})
xlabel('time','Fontsize',6)
ylabel('kicksize','Fontsize',6)
set(gca,'FontSize',6)
%link axes
linkaxes([ax ax2 ax3]);
set([ax2 ax3], 'XLim', xlim(ax));
set([ax2 ax3], 'Color', 'none');

%% save to the right size
set(gcf,'PaperUnits','inches');
oldsizes = get(gcf,'PaperPosition');
% This returns [x y width height]
newwidth = 3.2;
newheight = 1/2*oldsizes(4)/oldsizes(3)*newwidth;
set(gcf,'PaperPosition',[0 0 newwidth newheight]);
print('-opengl','6A.eps','-depsc','-r300')