close all;

% Sort the trials by frequency (needed for countorf function)
[f_star_sorted,sort_index] = sortrows(f_star_commanded);
U_star_sorted = 1./f_star_sorted;
A_star_sorted = A_star_measured(sort_index,:);
powercoef_mean_sorted = powercoef_mean(sort_index,:);
delay_sorted = delay(sort_index,:);
% powercoef_conv = powercoef_convtest(1,sort_index,:);

% [gradf,gradA] = gradient(powercoef_mean_sorted);

% Plot acceleration limit
acc_limit = 4.9; % Acceleration limit in m/s^2
v_limit = 0.5;

hold on

% Plot Cp vs. A* and f*
contourf(f_star_sorted,A_star_sorted,powercoef_mean_sorted,120,'LineStyle','none') %,[],'LineStyle','none'
% contourf(f_star_sorted,A_star_sorted,squeeze(powercoef_conv(1,:,:))./powercoef_mean_sorted,120,'LineStyle','none') %,[],'LineStyle','none'
% quiver(f_star_sorted,A_star_sorted,gradA,gradf)

caxis([-0.2 0.2])
% colorbarpwn(-6.0,0.2,'colorN',[0 0.5 1],'log',1.5)
colormap(bluewhitered)

contour(f_star_sorted,A_star_sorted,powercoef_mean_sorted,[0,0],'LineWidth',3,'LineColor','k','LineStyle','--')
scatter(f_star_sorted,A_star_sorted,[],'.','k')
grid on
xlabel('f* = fD/U')
ylabel('A* = A/D')
xlim([0.057 0.245])
ylim([-0.04 1.12])
set(gca, 'Layer', 'top')
grid off
% fstarvalues_extended = (0:0.025:0.3);
% a_limit_curve = acc_limit./(chord*(2*pi*(U.*fstarvalues_extended/chord)).^2);
% v_limit_curve = v_limit./(chord*(2*pi*(U.*fstarvalues_extended/chord)));
% plot(fstarvalues_extended,a_limit_curve)
% plot(fstarvalues_extended,v_limit_curve)

% % Plot Cp vs. A* and U*
% % contourf(U_star_sorted,A_star_sorted,powercoef_mean_sorted,120,'LineStyle','none') %,[],'LineStyle','none'
% % contourf(f_star_sorted,A_star_sorted,powercoef_mean_sorted,120,'LineStyle','none') %,[],'LineStyle','none'
% % quiver(f_star_sorted,A_star_sorted,gradA,gradf)
% caxis([-1.0 0.2])
% colormap(bluewhitered)
% contour(U_star_sorted,A_star_sorted,powercoef_mean_sorted,[0,0],'LineWidth',3,'LineColor','k','LineStyle','--')
% % scatter(U_star_sorted,A_star_sorted,[],'k')
% grid on
% xlabel('{\it U}* = {\it U/Df }')
% ylabel('{\it A}* = {\it A/D}')
% xlim([3 20])
% ylim([0 1.2])
% Ustarvalues_extended = (0:1:11);
% a_limit_curve = acc_limit./(chord*(2*pi*(flowspeed_fixed./(chord*Ustarvalues_extended))).^2);
% v_limit_curve = v_limit./(chord*(2*pi*(flowspeed_fixed./(chord*Ustarvalues_extended))));
% plot(Ustarvalues_extended,a_limit_curve)
% plot(Ustarvalues_extended,v_limit_curve)

c=colorbar();
c.Label.String = '$C_p$';
c.Label.Interpreter = 'Latex';



hold off