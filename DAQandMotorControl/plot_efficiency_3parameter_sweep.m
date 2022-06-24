%% Plotting uncorrected efficiency

clear;

load('20220619_TandemFoil_efficiency_A2E_a15_PHPh.mat');
% load('20220619_TandemFoil_efficiency_A2E_a33_PHPh.mat');
% load('20220619_TandemFoil_efficiency_A2E_a68_PHPh.mat');

lvlstp = 0.01;

[X, Y] = meshgrid(h3,p3);


%%
figure(1)

Eff_phys_2_mean = mean(Eff_phys_2,'all');
% Eff_phys_2_std = std(Eff_phys_2,'all');
Eff_phys_2_var = 100*abs((Eff_phys_2_mean-Eff_phys_2)/Eff_phys_2_mean);

sgtitle('Leading foil, percentage deviation from mean efficiency')
colormap('winter')

for i = 1:size(Eff_phys_2,2)
    
    Eff2 = squeeze(Eff_phys_2_var(:,i,:));

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

sgtitle('Uncorrected')
colormap('turbo')

for i = 1:size(Eff_phys_2,2)
    
    Eff3 = squeeze(Eff_phys_3(:,i,:));

    subplot(2,4,i);
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none');
    colorbar();
    caxis([min(Eff_phys_3,[],'all'), max(Eff_phys_3,[],'all')]);
    
    t = ['phase = ',num2str(ph(i))];
    title(t);
    
end

%% Corrected efficiency

% figure(3)
% 
% sgtitle('Corrected')
% colormap('turbo')
% 
% for i = 1:size(Eff_corr_2,2)
%     
%     Eff2 = squeeze(Eff_corr_2(:,i,:));
% 
%     subplot(2,4,i);
%     contourf(X, Y, Eff2, 'LevelStep', lvlstp, 'LineStyle', 'none');
%     colorbar();
%     caxis([0.23, 0.246]);
%     
%     t = ['phase = ',num2str(ph(i))];
%     title(t);
%     
% end
% 
% %%
% figure(4)
% 
% sgtitle('Corrected')
% colormap('turbo')
% 
% for i = 1:size(Eff_corr_2,2)
%     
%     Eff3 = squeeze(Eff_corr_3(:,i,:));
% 
%     subplot(2,4,i);
%     contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none');
%     colorbar();
%     caxis([0.05, 0.22]);
%     
%     t = ['phase = ',num2str(ph(i))];
%     title(t);
%     
% end