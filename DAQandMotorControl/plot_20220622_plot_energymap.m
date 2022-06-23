figure(1)

h3 = [0.5500    0.7000    0.8500    1.0000    1.1500    1.3000];
p3 = [60    65    70    75    80];
ph = [-180  -120   -60     0    60   120   180];

[X, Y] = meshgrid(p3, h3);
lvlstp = 0.005;

for i = 1:7
    subplot(2,4,i);
    Eff3 = squeeze(Eff_3(i,:,:));
    
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none');
    colorbar();
    
    xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel('$h_0$','Interpreter','latex', 'FontSize', 20);
    set(gca,'FontSize',16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
    t = ['phase = ', num2str(ph(i))];
    title(t, 'FontSize',16, 'LineWidth', 1.5);
    caxis([-0.03, 0.21]);
    
    colormap('turbo')
    
end
