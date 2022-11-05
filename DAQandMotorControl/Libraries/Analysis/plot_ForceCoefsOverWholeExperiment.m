close all
hold on;
plot(reshape(liftcoef_alltrials',[1,numel(liftcoef_alltrials)]),'DisplayName','{\it C}_L','Color','red','LineWidth',4)
plot(reshape(dragcoef_alltrials',[1,numel(dragcoef_alltrials)]),'DisplayName','{\it C}_D','Color','blue','LineWidth',4)
% plot(liftcoef_alltrials(:,1),'DisplayName','{\it C}_L','Color','red','LineWidth',4)
% plot(dragcoef_alltrials(:,1),'DisplayName','{\it C}_D','Color','blue','LineWidth',4)
ylim([-1 4])
% xlim([0 17])
legend()
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2);
xlabel('Trial number (total time 13.3 hours)')
hold off;