%% Generating a single foil gif

% working directory:
cd('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=68_p3=75_h3=0.8_ph=60_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\tandem_foil');

% load compiled data file:
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=68_p3=75_h3=0.8_ph=60_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.68_p3=75_h3=0.8c_ph=60.mat');

frame = 23;

chord = 0.061;
U = 0.33;

for n = 1:frame
    
    h = figure;
    set(h,'color','w','unit','centimeters','position',[2,2,40,20]); % for vorticity and forces
    
%     frame = 2*(n + 2);
    x = compiled_data(n).x+3;
    y = compiled_data(n).y;
    foil2_coords = compiled_data(n).foil2_coords;
    foil3_coords = compiled_data(n).foil3_coords;
    vort = compiled_data(n).vort*chord/U;
    timeStep = compiled_data(n).timeStep;
    
    isValid = compiled_data(n).isValid;
    vort(isValid == 0) = NaN;
    
    contourf(x,y,vort,'LineStyle','none','LevelStep',0.1); hold on;
    axis equal
    caxis([-2,2])
    ylim([-2,2])
%     xlim([-1,1])
    
    plot(foil2_coords([1,3,5])+3, foil2_coords([2,4,6]), 'k', 'LineWidth', 10); % plot leading foil
    plot(foil3_coords([1,3,5])+3, foil3_coords([2,4,6]), 'k', 'LineWidth', 10); hold off; % plot trailing foil
    
    xlabel('$x/c$','fontsize',20,'interpreter','Latex');
    ylabel('$y/c$','fontsize',20,'interpreter','Latex');
    set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in');
    
    colormap(brewermap([],"-RdBu"))
    
    clrbr = colorbar('Linewidth',1.5);
    clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
    clrbr.Label.Interpreter = 'latex';
    
    figname = ['fig', num2str(n), '.tif'];
%     saveas(gcf,figname) % alternate way of saving
    export_fig(figname); % saves the .tif file
    
end

%% Generating the gif

gifname = '20221006_TandemFoils_aT4=0.68_p3=75_h3=0.8c_ph=60.gif';

freq_piv = frame*(0.12*U/chord);

for ii = 1:frame % 23 is the number of frames per cycle (14... Hz/0.6.. Hz)
    figname = ['fig' num2str(ii) '.tif'];
    Z = imread(figname);
    [imind,cm] = rgb2ind(Z,256); 
    
    if ii == 1 
        imwrite(imind,cm,gifname,'gif','Loopcount',inf,'DelayTime',1/freq_piv); 
    else 
        imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',1/freq_piv); 
    end
end

