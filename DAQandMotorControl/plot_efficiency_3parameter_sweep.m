%% Plotting uncorrected efficiency

clear;

load('20220619_TandemFoil_efficiency_A2E_a15_PHPh.mat');
% load('20220619_TandemFoil_efficiency_A2E_a33_PHPh.mat');
% load('20220619_TandemFoil_efficiency_A2E_a68_PHPh.mat');

lvlstp = 0.01;

[X, Y] = meshgrid(p3,h3);

%%
figure(1)

Eff_phys_2_mean = mean(Eff_phys_2,'all');
Eff_phys_2_var = 100*abs((Eff_phys_2_mean-Eff_phys_2)/Eff_phys_2_mean);

sgtitle('Leading foil, percentage deviation from mean efficiency')
colormap('viridis')

for i = 1:size(Eff_phys_2,2)
    
    Eff2 = squeeze(Eff_phys_2_var(:,i,:))';

    subplot(2,4,i);
    contourf(X, Y, Eff2, 'LineStyle', 'none');%, 'LevelStep', 0.01);
    c = colorbar();
	caxis([0, max(Eff_phys_2_var,[],'all')]);
    c.Label.String = ('|100*(eff_mean-eff)/eff_mean|');
    
    t = ['phase = ',num2str(ph(i))];
    title(t);
    
end

%%
figure(2)

sgtitle('Uncorrected, trailing', 'FontSize', 20, 'Interpreter', 'latex');
colormap('magma')

for i = 1:size(Eff_phys_2,2)
    
    Eff3 = squeeze(Eff_phys_3(:,i,:))';

    subplot(3,3,i);
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;
    
    % Find maximum
    
    [M, I] = max(Eff3,[],'all','linear');
%     xline(X(I),'--y');
%     yline(Y(I),'--y');
    plot(X(I),Y(I),'p','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','y');
    
    
    colorbar();
    caxis([min(Eff_phys_3,[],'all'), max(Eff_phys_3,[],'all')]);
    
    t = ['phase = ',num2str(ph(i)), '$^o$, $\eta_{tr}$ = ', num2str(M,3)];
    title(t, 'Interpreter', 'latex', 'FontSize', 18);
    
end

%% Corrected efficiency
% 
% figure(3)
% 
% sgtitle('Corrected')
% colormap('turbo')
% 
% for i = 1:size(Eff_corr_2,2)
%     
%     Eff2 = squeeze(Eff_corr_2(:,i,:))';
% 
%     subplot(2,4,i);
%     contourf(X, Y, Eff2, 'LevelStep', lvlstp, 'LineStyle', 'none');
%     colorbar();
% %     caxis([0.23, 0.246]);
%     
%     t = ['phase = ',num2str(ph(i))];
%     title(t);
%     
% end
% 
%%
figure(4)

sgtitle('Corrected, trailing', 'FontSize', 20, 'Interpreter', 'latex');
colormap('magma')

for i = 1:size(Eff_corr_2,2)
    
    Eff3 = squeeze(Eff_corr_3(:,i,:))';

    subplot(3,3,i);
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;
    
    % Find maximum
    
    [M, I] = max(Eff3,[],'all','linear');
%     xline(X(I),'--y');
%     yline(Y(I),'--y');
    plot(X(I),Y(I),'p','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','y');
    
    
    colorbar();
    caxis([0, max(Eff_corr_3,[],'all')]);
    
    t = ['phase = ',num2str(ph(i)), '$^o$, $\eta_{tr}$ = ', num2str(M,3)];
    title(t, 'Interpreter', 'latex', 'FontSize', 18);
    
end

%%
figure(5)

sgtitle('System Efficiency', 'FontSize', 20, 'Interpreter', 'latex');
colormap('magma')

for i = 1:size(Eff_sys,2)
    
    Eff3 = squeeze(Eff_sys(:,i,:))';

    subplot(3,3,i);
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;
    
    % Find maximum
    
    [M, I] = max(Eff3,[],'all','linear');
%     xline(X(I),'--y');
%     yline(Y(I),'--y');
    plot(X(I),Y(I),'p','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','y');
    
    
    colorbar();
    caxis([min(Eff_sys,[],'all'), max(Eff_sys,[],'all')]);
    
    t = ['phase = ',num2str(ph(i)), '$^o$, $\eta_{sys}$ = ', num2str(M,3)];
    title(t, 'Interpreter', 'latex', 'FontSize', 18);
    
end