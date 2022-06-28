%% Plotting uncorrected efficiency

clear;

% load('20220619_TandemFoil_efficiency_A2E_a15_PHPh.mat');
load('20220619_TandemFoil_efficiency_A2E_a33_PHPh.mat');
% load('20220619_TandemFoil_efficiency_A2E_a68_PHPh.mat');

lvlstp = 0.01;

[X, Y] = meshgrid(p3,h3);

%%
figure(1)

colormap('viridis')

Eff_phys_2_std = std(Eff_phys_2,0,2);
Eff_phys_2_mean = mean(Eff_phys_2,'all');
Eff2 = squeeze(Eff_phys_2_std(:,1,:))';

contourf(X, Y, Eff2, 'LineStyle', 'none');

xlabel('$\theta_{tr}$', 'Interpreter', 'latex');
ylabel('$h_{tr}$', 'Interpreter', 'latex');

c = colorbar();
caxis([0, max(Eff_phys_2_std,[],'all')]);
set(gca,'FontSize',22,'TickLabelInterpreter', 'latex');

t = ['STD($\eta_{le}$), $\bar \eta_{le} =$ ', num2str((Eff_phys_2_mean),3)];
title(t, 'Interpreter', 'latex', 'FontSize', 24);

%%
figure(2)

mt = ['Unorrected trailing foil efficiency, $\alpha_{T/4}$ = ', num2str(alphaT4,2)];
sgtitle(mt, 'FontSize', 20, 'Interpreter', 'latex');
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

mt = ['Corrected trailing foil efficiency, $\alpha_{T/4}$ = ', num2str(alphaT4,2)];
sgtitle(mt, 'FontSize', 20, 'Interpreter', 'latex');
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

mt = ['System efficiency, $\alpha_{T/4}$ = ', num2str(alphaT4,2)];
sgtitle(mt, 'FontSize', 20, 'Interpreter', 'latex');
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