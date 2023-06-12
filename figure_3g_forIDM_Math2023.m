%% Figure 3g - no excursion after a single kick
% Generate the separation curve

paramset0 = [0.35,20,0.2,0.04,5e-7,2e-5,0.8,7e7,3,1*10^-7,2,1e-4,0.01,12000,0.2,10,3];
fig3 = figure_3g_onekickexcursioncalculationMay2023;

IC = zeros(6,1);
IC(2) = 7e7;
IC(1) = 4000;

% Start all runs with an initial kick of 4000 and a 2nd kick of a different
% size.

newtimes = [5:10,11:2:30,31:5:81,90:10:410];
kickestabove = newtimes*0+270000; % above the threshold at all times

[timelist,kicklist] = fig3.curve(newtimes,IC,kickestabove,1,paramset0);
kb(:,1) = kicklist(:,1);
ka(:,1) = kicklist(:,2);

%% Figure 3g - plot to confirm excursions / no excursions in the data
% Check four at a time.  Pause after each set of four.
for k3 = 1:floor(length(timelist)/4)
    for k1 = 1:4
        k2 = k1+2+4*(k3-1);
        subplot(2,4,k1)
        fig3.plot(IC,ka(k2),timelist(k2),paramset0)
        title(['excursion, tau = ',num2str(timelist(k2))])
        subplot(2,4,k1+4)
        fig3.plot(IC,kb(k2),timelist(k2),paramset0)
         title(['no excursion, tau = ',num2str(timelist(k2))])
    end
    pause
end

%% Figure 3g - create part G of the figure
figure(1)

scale = 1e5;
ax = gca;
% Use the "below" curve for the separation curve.
area(timelist,kb(:,1)/scale)
hold on
% Mark locations used in the time series plots
%plot(283,15000/scale,'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',2)
%plot(283,14800/scale,'o','MarkerFaceColor',"#EDB120",'MarkerEdgeColor','k','MarkerSize',2)
axis([0 400 0 3])
xticks([0:50:400])
ax.XAxis.FontSize = 6;
ax.YAxis.FontSize = 6;
yticks([0:.5:3])
%yticklabels({'0','0.2','0.4','0.6','0.8','1',''})
grid on
xlabel("exposure period",'Fontsize',8)
ylabel("viral exposure x1e5",'Fontsize',8)
set(gca,'LooseInset',get(gca,'TightInset'));
text(375,2.75,'G','EdgeColor','k','fontsize',8,'Margin',0.5)
hold off

eps_save(1,'fig3g-secondkick.eps')


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