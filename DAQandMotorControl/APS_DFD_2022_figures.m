%% APS DFD plots

clear;

%% 5. Wake regimes - three aT4

chord = 0.061;
U = 0.33;

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=16_p3=75_h3=0.7_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU_01\Export\TandemFoil_aT4=0.16_p3=75_h3=0.7c_ph=-120.mat');
n = 11;
x_low = compiled_data(n).x+3;
y_low = compiled_data(n).y;
vort_low = compiled_data(n).vort*chord/U;
isValid = compiled_data(n).isValid;
vort_low(isValid == 0) = NaN;

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=33_p3=75_h3=0.7_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.33_p3=75_h3=0.7c_ph=-120.mat');
n = 11;
x_mid = compiled_data(n).x+3;
y_mid = compiled_data(n).y;
vort_mid = compiled_data(n).vort*chord/U;
isValid = compiled_data(n).isValid;
vort_mid(isValid == 0) = NaN;

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=68_p3=75_h3=0.7_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.68_p3=75_h3=0.7c_ph=-120.mat');
n = 23;
x_hig = compiled_data(n).x+3;
y_hig = compiled_data(n).y;
vort_hig = compiled_data(n).vort*chord/U;
isValid = compiled_data(n).isValid;
vort_hig(isValid == 0) = NaN;


h = figure(1);
set(h,'color','w');

subplot(3,1,1)
% title('$\alpha_{T/4} = 0.16$','interpreter','latex');
contourf(x_low,y_low,vort_low,'LineStyle','none','LevelStep',0.1);
colormap(brewermap([],"-RdBu"));
axis equal
caxis([-1,1]);
xlim([-0.8,5]);
ylim([-2,2]);

ylabel('$y/c$','fontsize',22,'interpreter','Latex');
set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','XTick',[]);

% clrbr = colorbar('Linewidth',1.5);
% clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
% clrbr.Label.Interpreter = 'latex';


subplot(3,1,2)
% title('$\alpha_{T/4} = 0.33$','interpreter','latex');
contourf(x_mid,y_mid,vort_mid,'LineStyle','none','LevelStep',0.05);
colormap(brewermap([],"-RdBu"));
axis equal
caxis([-1,1]);
xlim([-0.8,5]);
ylim([-2,2]);

ylabel('$y/c$','fontsize',22,'interpreter','Latex');
set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','XTick',[]);

% clrbr = colorbar('Linewidth',1.5);
% clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
% clrbr.Label.Interpreter = 'latex';


subplot(3,1,3)
% title('$\alpha_{T/4} = 0.68$','interpreter','latex');
contourf(x_hig,y_hig,vort_hig,'LineStyle','none','LevelStep',0.1);
colormap(brewermap([],"-RdBu"));
axis equal
caxis([-1,1]);
xlim([-0.8,5]);
ylim([-2,2]);

ylabel('$y/c$','fontsize',22,'interpreter','Latex');
xlabel('$x/c$','fontsize',22,'interpreter','Latex');
set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in');

clrbr = colorbar('Linewidth',1.5);
clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
clrbr.Label.Interpreter = 'latex';

%% 6. Tandem Foils

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=33_p3=75_h3=0.7_ph=60_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.33_p3=75_h3=0.7c_ph=60.mat');

chord = 0.061;
U = 0.33;

h = figure(2);
set(h,'color','w');

tshot = [0.0,0.3,0.6];
frameS = [17, 1, 8, 13];

for n = 1:3
    
%     frame = 11 + 4*n;
    frame = frameS(n);
    x = compiled_data(frame).x+3;
    y = compiled_data(frame).y;
    foil2_coords = compiled_data(frame).foil2_coords;
    foil3_coords = compiled_data(frame).foil3_coords;
    vort = compiled_data(frame).vort*chord/U;
    timeStep = compiled_data(frame).timeStep;
    
    isValid = compiled_data(frame).isValid;
    vort(isValid == 0) = NaN;
    
    subplot(1,3,n)
    contourf(x,y,vort,'LineStyle','none','LevelStep',0.1); hold on;
    axis equal
    caxis([-1,1])
    
    plot(foil2_coords([1,3,5])+3, foil2_coords([2,4,6]), 'k', 'LineWidth', 5); % plot leading foil
    plot(foil3_coords([1,3,5])+3, foil3_coords([2,4,6]), 'k', 'LineWidth', 5); hold off; % plot trailing foil
    
    xlabel('$x/c$','fontsize',16,'interpreter','Latex');
    if n == 1
        ylabel('$y/c$','fontsize',16,'interpreter','Latex');
        set(gca,'fontname','Times New Roman','fontsize',16,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in');
    else
        set(gca,'fontname','Times New Roman','fontsize',16,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','YTick',[]);
    end
%     title(['t/T = ',num2str(round(timeStep*(1/1000)/(1/0.649),1))],'Interpreter','latex');
    title(['t/T = ',num2str(tshot(n))],'Interpreter','latex');
    
end

clrbr = colorbar('Linewidth',1.5);
clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
clrbr.Label.Interpreter = 'latex';
colormap(brewermap([],"-RdBu"))

%% 7. pitch for all aT4 cases

clear;

% single
load('20221011_SingleFoil_efficiency_A3E_PH.mat');
Eff_2_single = Eff_2;
Eff_2_std_single = Eff_2_std;
% tandem
load('20221011_TandemFoil_efficiency_A3E_a155_330_679_PHPh_CpFrstrm_EffFrstrm_SysEffFrstrm_wBaseline.mat');

aT4 = [0.155, 0.33, 0.679]; laT4 = length(aT4);
p3 = [65  70  75]; lp3 = length(p3);
h3 = [0.7000    0.8000    0.9000    1.0000    1.1000    1.2000]; lh3 = length(h3);
ph = [-180  -120   -60     0    60   120]; lph = length(ph);

for n = 1:lp3
    eff_2_single(n) = Eff_2_single(1,n);
    eff_2_std_single(n) = Eff_2_std_single(1,n);
    for m = 1:laT4
        eff_3(m,n) = Eff_3(m,n,1,2);
        eff_3_std(m,n) = Eff_3_std(m,n,1,2);
        eff_sys(m,n) = Eff_sys(m,n,1,2);
        eff_sys_std(m,n) = Eff_sys_std(m,n,1,2);
    end
end

figure(3)

subplot(1,3,1)
errorbar(p3, eff_2_single, eff_2_std_single,'-','linewidth',1.5,'color',[0.2,0,1]);
title('Single Foil','interpreter','latex');
xlabel('$\theta_{0}$','interpreter','latex');
ylabel('$\eta$','interpreter','latex');
legend('Single foil, $h_{0} = 0.7c$','interpreter','latex','location','south');
ylim([0,0.6]);
xlim([60,80]);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');

subplot(1,3,2)
errorbar(p3, eff_3(1,:), eff_3_std(1,:),'-','linewidth',1.5); hold on;
errorbar(p3, eff_3(2,:), eff_3_std(2,:),'-','linewidth',1.5);
errorbar(p3, eff_3(3,:), eff_3_std(3,:),'-','linewidth',1.5,'color',[0.4660 0.6740 0.1880]); hold off;
title('Trailing Foil','interpreter','latex');
xlabel('$\theta_{0,\mathrm{tr}}$','interpreter','latex');
ylabel('$\eta_{\mathrm{tr}}$','interpreter','latex');
ylim([0,0.6]);
xlim([60,80]);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');

subplot(1,3,3)
errorbar(p3, eff_sys(1,:), eff_sys_std(1,:),'-','linewidth',1.5); hold on;
errorbar(p3, eff_sys(2,:), eff_sys_std(2,:),'-','linewidth',1.5);
errorbar(p3, eff_sys(3,:), eff_sys_std(3,:),'-','linewidth',1.5,'color',[0.4660 0.6740 0.1880]); hold off;
title('System','interpreter','latex');
xlabel('$\theta_{0,\mathrm{tr}}$','interpreter','latex');
ylabel('$\eta_{\mathrm{sys}}$','interpreter','latex');
legend('$\alpha_{T/4} = 0.16$, $h_{0,\mathrm{tr}} = 0.7c$',...
       '$\alpha_{T/4} = 0.33$, $h_{0,\mathrm{tr}} = 0.7c$',...
       '$\alpha_{T/4} = 0.68$, $h_{0,\mathrm{tr}} = 0.7c$',...
    'interpreter','latex','location','south');
ylim([0,0.6]);
xlim([60,80]);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');

%% 8. Heave for all aT4 cases

clear;

% single
load('20221011_SingleFoil_efficiency_A3E_PH.mat');
Eff_2_single = Eff_2;
Eff_2_std_single = Eff_2_std;
% tandem
load('20221011_TandemFoil_efficiency_A3E_a155_330_679_PHPh_CpFrstrm_EffFrstrm_SysEffFrstrm_wBaseline.mat');

aT4 = [0.155, 0.33, 0.679]; laT4 = length(aT4);
p3 = [65  70  75]; lp3 = length(p3);
h3 = [0.7000    0.8000    0.9000    1.0000    1.1000    1.2000]; lh3 = length(h3);
ph = [-180  -120   -60     0    60   120]; lph = length(ph);

for n = 1:lh3
    eff_2_single(n) = Eff_2_single(n,3);
    eff_2_std_single(n) = Eff_2_std_single(n,3);
    for m = 1:laT4
        eff_3(m,n) = Eff_3(m,3,n,2);
        eff_3_std(m,n) = Eff_3_std(m,3,n,2);
        for q = 1:lph
            eff_sys(m,n,q) = Eff_sys(m,3,n,q);
            eff_sys_std(m,n,q) = Eff_sys_std(m,3,n,q);
        end
    end
end

figure(4)

subplot(1,3,1)
errorbar(h3, eff_2_single, eff_2_std_single,'-','linewidth',1.5,'color',[0.2,0,1]);
title('Single Foil','interpreter','latex');
xlabel('$h_{0}/c$','interpreter','latex');
ylabel('$\eta$','interpreter','latex');
legend('Single foil, $\theta_{0} = 75^{\circ}$','interpreter','latex','location','south');
ylim([0,0.6]);
xlim([0.6,1.3]);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');

subplot(1,3,2)
errorbar(h3, eff_3(1,:), eff_3_std(1,:),'-','linewidth',1.5); hold on;
errorbar(h3, eff_3(2,:), eff_3_std(2,:),'-','linewidth',1.5);
errorbar(h3, eff_3(3,:), eff_3_std(3,:),'-','linewidth',1.5,'color',[0.4660 0.6740 0.1880]); hold off;
title('Trailing Foil','interpreter','latex');
xlabel('$h_{0,\mathrm{tr}}/c$','interpreter','latex');
ylabel('$\eta_{\mathrm{tr}}$','interpreter','latex');
ylim([0,0.6]);
xlim([0.6,1.3]);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');

subplot(1,3,3)
errorbar(h3, eff_sys(1,:,2), eff_sys_std(1,:,2),'-','linewidth',1.5); hold on;
errorbar(h3, eff_sys(2,:,2), eff_sys_std(2,:,2),'-','linewidth',1.5);
errorbar(h3, eff_sys(3,:,2), eff_sys_std(3,:,2),'-','linewidth',1.5,'color',[0.4660 0.6740 0.1880]); hold off;
title('System','interpreter','latex');
xlabel('$h_{0,\mathrm{tr}}/c$','interpreter','latex');
ylabel('$\eta_{\mathrm{sys}}$','interpreter','latex');
legend('$\alpha_{T/4} = 0.16$, $\theta_{0,\mathrm{tr}} = 75^{\circ}$',...
       '$\alpha_{T/4} = 0.33$, $\theta_{0,\mathrm{tr}} = 75^{\circ}$',...
       '$\alpha_{T/4} = 0.68$, $\theta_{0,\mathrm{tr}} = 75^{\circ}$',...
    'interpreter','latex','location','south');
ylim([0,0.6]);
xlim([0.6,1.3]);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');

figure(5)

errorbar(h3, eff_sys(3,:,1), eff_sys_std(3,:,1),'-','linewidth',1.5); hold on;
errorbar(h3, eff_sys(3,:,2), eff_sys_std(3,:,2),'-','linewidth',1.5);
errorbar(h3, eff_sys(3,:,3), eff_sys_std(3,:,3),'-','linewidth',1.5);
errorbar(h3, eff_sys(3,:,4), eff_sys_std(3,:,4),'-','linewidth',1.5);
errorbar(h3, eff_sys(3,:,5), eff_sys_std(3,:,5),'-','linewidth',1.5);
errorbar(h3, eff_sys(3,:,6), eff_sys_std(3,:,6),'-','linewidth',1.5); hold off;
title('System','interpreter','latex');
xlabel('$h_{0,\mathrm{tr}}/c$','interpreter','latex');
ylabel('$\eta_{\mathrm{sys}}$','interpreter','latex');
% legend('$\alpha_{T/4} = 0.16$, $\theta_{0,\mathrm{tr}} = 75^{\circ}$',...
%        '$\alpha_{T/4} = 0.33$, $\theta_{0,\mathrm{tr}} = 75^{\circ}$',...
%        '$\alpha_{T/4} = 0.68$, $\theta_{0,\mathrm{tr}} = 75^{\circ}$',...
%     'interpreter','latex','location','south');
ylim([0,0.6]);
xlim([0.6,1.3]);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');

%% 9. Phase for all aT4 cases

clear;

% single
load('20221011_SingleFoil_efficiency_A3E_PH.mat');
Eff_2_single = Eff_2;
Eff_2_std_single = Eff_2_std;
% tandem
load('20221011_TandemFoil_efficiency_A3E_a155_330_679_PHPh_CpFrstrm_EffFrstrm_SysEffFrstrm_wBaseline.mat');

Eff_3(:,:,:,7) = Eff_3(:,:,:,1); % repeat the -180 datapoint for the +180 datapoint
Eff_3_std(:,:,:,7) = Eff_3_std(:,:,:,1);
Eff_sys(:,:,:,7) = Eff_sys(:,:,:,1);
Eff_sys_std(:,:,:,7) = Eff_sys_std(:,:,:,1);

aT4 = [0.155, 0.33, 0.679]; laT4 = length(aT4);
p3 = [65  70  75]; lp3 = length(p3);
h3 = [0.7000    0.8000    0.9000    1.0000    1.1000    1.2000]; lh3 = length(h3);
ph = [-180  -120   -60     0    60   120   180]; lph = length(ph); % repeating the last phase datapoints

for q = 1:lh3
    for n = 1:lph
        for m = 1:laT4
            eff_3(m,n,q) = Eff_3(m,3,q,n);
            eff_3_std(m,n,q) = Eff_3_std(m,3,q,n);
            eff_sys(m,n,q) = Eff_sys(m,3,q,n);
            eff_sys_std(m,n,q) = Eff_sys_std(m,3,q,n);
        end % (aT4, PH, H3)
    end
end

figure(4)

subplot(1,2,1)
errorbar(ph, eff_3(1,:,1), eff_3_std(1,:,1),'-','linewidth',1.5); hold on;
errorbar(ph, eff_3(2,:,1), eff_3_std(2,:,1),'-','linewidth',1.5);
errorbar(ph, eff_3(3,:,1), eff_3_std(3,:,1),'-','linewidth',1.5,'color',[0.4660 0.6740 0.1880]); hold off;
title('Trailing Foil','interpreter','latex');
xlabel('$\psi_{1-2}$','interpreter','latex');
ylabel('$\eta_{\mathrm{tr}}$','interpreter','latex');
ylim([0,0.6]);
xlim([-180,180]);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');

subplot(1,2,2)
errorbar(ph, eff_sys(1,:,1), eff_sys_std(1,:,1),'-','linewidth',1.5); hold on;
errorbar(ph, eff_sys(2,:,1), eff_sys_std(2,:,1),'-','linewidth',1.5);
errorbar(ph, eff_sys(3,:,1), eff_sys_std(3,:,1),'-','linewidth',1.5,'color',[0.4660 0.6740 0.1880]); hold off;
title('System','interpreter','latex');
xlabel('$\psi_{1-2}$','interpreter','latex');
ylabel('$\eta_{\mathrm{sys}}$','interpreter','latex');
legend('$\alpha_{T/4} = 0.16$, $\theta_{0,\mathrm{tr}} = 75^{\circ}$, $h_{0,\mathrm{tr}} = 0.7c$',...
       '$\alpha_{T/4} = 0.33$, $\theta_{0,\mathrm{tr}} = 75^{\circ}$, $h_{0,\mathrm{tr}} = 0.7c$',...
       '$\alpha_{T/4} = 0.68$, $\theta_{0,\mathrm{tr}} = 75^{\circ}$, $h_{0,\mathrm{tr}} = 0.7c$',...
    'interpreter','latex','location','south');
ylim([0,0.6]);
xlim([-180,180]);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');


figure(5)
hold on;
for q = 1:lh3
    % errorbar(ph, eff_3(1,:,1), eff_3_std(1,:,1),'-','linewidth',1.5); hold on;
    % errorbar(ph, eff_3(2,:,1), eff_3_std(2,:,1),'-','linewidth',1.5);
    plot(ph, eff_3(3,:,q),'-','linewidth',1.5); %hold off;
end
hold off;
title('Trailing Foil','interpreter','latex');
xlabel('$\psi_{1-2}$','interpreter','latex');
ylabel('$\eta_{\mathrm{tr}}$','interpreter','latex');
legend('h = 0.7c','h = 0.8c','h = 0.9c','h = 1.0c','h = 1.1c','h = 1.2c');
ylim([0,0.2]);
xlim([-180,180]);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');

%% 10. LEV trajectory

clear;

frame = 23;
half = 11;

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=33_p3=75_h3=0.7_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.33_p3=75_h3=0.7c_ph=-120.mat');

i = 0;

for n = [12:23,1:3]
    i = i + 1;
    x_33(:,:,i) = compiled_data(n).x;
    y_33(:,:,i) = compiled_data(n).y;
    foil2_coords_33(:,i) = compiled_data(n).foil2_coords;
    foil3_coords_33(:,i) = compiled_data(n).foil3_coords;
    u_33(:,:,i) = compiled_data(n).u;
    v_33(:,:,i) = compiled_data(n).v;
    vort = compiled_data(n).vort;
    isValid = compiled_data(n).isValid;
    vort(isValid == 0) = NaN; % foil mask
    vort_33(:,:,i) = vort;
    timeIndex_33(i) = compiled_data(n).timeStep + 1;
    
    ds = x_33(1,1,i) - x_33(1,2,i); % distance between each vector
    [dudx, dudy] = gradient(u_33(:,:,i), ds);
    [dvdx, dvdy] = gradient(v_33(:,:,i), ds);

    Q_mid = (dudx.*dvdy - dvdx.*dudy);
    Q_mid_neg = zeros(size(Q_mid));
    Q_mid_neg(Q_mid>1) = Q_mid(Q_mid>1);
    Q_mid_neg(log(Q_mid)<13) = 0;

    Q_mid_neg(isValid == 0) = NaN; % foil mask

    Q_mid = log(Q_mid_neg);
    Q_mid(Q_mid == -Inf) = 0;
    Q_33(:,:,i) = Q_mid;
    
end

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=68_p3=75_h3=0.7_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.68_p3=75_h3=0.7c_ph=-120.mat');

i = 0;

for n = [1:15]
    i = i + 1;
    x_68(:,:,i) = compiled_data(n).x;
    y_68(:,:,i) = compiled_data(n).y;
    foil2_coords_68(:,i) = compiled_data(n).foil2_coords;
    foil3_coords_68(:,i) = compiled_data(n).foil3_coords;
    u_68(:,:,i) = compiled_data(n).u;
    v_68(:,:,i) = compiled_data(n).v;
    vort = compiled_data(n).vort;
    isValid = compiled_data(n).isValid;
    vort(isValid == 0) = NaN; % foil mask
    vort_68(:,:,i) = vort;
    timeIndex_68(i) = compiled_data(n).timeStep + 1;
    
    ds = x_68(1,1,i) - x_68(1,2,i); % distance between each vector
    [dudx, dudy] = gradient(u_68(:,:,i), ds);
    [dvdx, dvdy] = gradient(v_68(:,:,i), ds);

    Q_mid = (dudx.*dvdy - dvdx.*dudy);
    Q_mid_neg = zeros(size(Q_mid));
    Q_mid_neg(Q_mid>1) = Q_mid(Q_mid>1);
    Q_mid_neg(log(Q_mid)<13) = 0;

    Q_mid_neg(isValid == 0) = NaN; % foil mask

    Q_mid = log(Q_mid_neg);
    Q_mid(Q_mid == -Inf) = 0;
    Q_68(:,:,i) = Q_mid;
    
end

% Plot trajectories and select manually vortex center

for do_it_again = 1:5
    for m = 1:15

        figure('name', 'a33', 'WindowState', 'maximized');
        contourf(x_33(:,:,m), y_33(:,:,m), Q_33(:,:,m),'LineStyle','none');
        axis equal
        caxis([0,20])
        xlim([-3.5,2.5])
        ylim([-2,2])
    %     colorbarpwn(0,20,'log','colorN',[1,1,1],'colorW',[0.2,0,0.8],'colorP',[1,0,1],'full', 5);
        LEV_ctr_33 = drawpoint();
        LEV_coords_33(do_it_again,m,:) = LEV_ctr_33.Position;
        close 'a33';

        figure('name', 'a68', 'WindowState', 'maximized');
        contourf(x_68(:,:,m), y_68(:,:,m), Q_68(:,:,m),'LineStyle','none');
        axis equal
        caxis([0,20])
        xlim([-3.5,2.5])
        ylim([-2,2])
    %     colorbarpwn(0,20,'log','colorN',[1,1,1],'colorW',[0.2,0,0.8],'colorP',[1,0,1],'full', 5);
        LEV_ctr_68 = drawpoint();
        LEV_coords_68(do_it_again,m,:) = LEV_ctr_68.Position;
        close 'a68';

    end
end

LEV_coords_33_avg = squeeze(mean(LEV_coords_33,1));
LEV_coords_68_avg = squeeze(mean(LEV_coords_68,1));

save('20221026_Main_vortex_tracking_a=0.33_p3=75_h3=0.7c_ph=-120.mat','LEV_coords_33','LEV_coords_33_avg');
save('20221026_Main_vortex_tracking_a=0.68_p3=75_h3=0.7c_ph=-120.mat','LEV_coords_68','LEV_coords_68_avg');

load('20221026_Main_vortex_tracking_a=0.33_p3=75_h3=0.7c_ph=-120.mat');
load('20221026_Main_vortex_tracking_a=0.68_p3=75_h3=0.7c_ph=-120.mat');

LEV_coords_33_avg(:,1) = LEV_coords_33_avg(:,1) + 3;
LEV_coords_68_avg(:,1) = LEV_coords_68_avg(:,1) + 3;

figure(7)

toime = 1:half;
plot(LEV_coords_33_avg(:,1),LEV_coords_33_avg(:,2),'^','color',[0.6350 0.0780 0.1840],'MarkerSize',10,'LineWidth',2,'MarkerFaceColor',[0.6350 0.0780 0.1840]); hold on;
plot(LEV_coords_68_avg(:,1),LEV_coords_68_avg(:,2),'^','color',[0.4660 0.6740 0.1880],'MarkerSize',10,'LineWidth',2,'MarkerFaceColor',[0.4660 0.6740 0.1880]); hold off;
legend('$\alpha_{T/4} = 0.33$','$\alpha_{T/4} = 0.68$','interpreter','latex','fontsize',22,'location','southeast');
tit = 'Primary Vortex Trajectory';
title(tit,'interpreter','latex','fontsize',24);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');
% axis equal
%%

figure(6)

subplot(1,2,1);
contourf(x_mid,y_mid,Q_mid,'LineStyle','none'); hold on;
axis equal
plot(foil2_coords_mid([1,3,5]), foil2_coords_mid([2,4,6]), 'k', 'LineWidth', 7); % plot leading foil
plot(foil3_coords_mid([1,3,5]), foil3_coords_mid([2,4,6]), 'k', 'LineWidth', 7); hold off; % plot trailing foil
caxis([0,20])
xlim([-3.5,2.5])
ylim([-2,2])
colorbarpwn(0,20,'log','colorN',[1,1,1],'colorW',[0.2,0,0.8],'colorP',[1,0,1],'full', 5);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in');

subplot(1,2,2);
contourf(x_hig,y_hig,Q_hig,'LineStyle','none'); hold on;
axis equal
plot(foil2_coords_hig([1,3,5]), foil2_coords_hig([2,4,6]), 'k', 'LineWidth', 7); % plot leading foil
plot(foil3_coords_hig([1,3,5]), foil3_coords_hig([2,4,6]), 'k', 'LineWidth', 7); hold off; % plot trailing foil
caxis([0,20])
xlim([-3.5,2.5])
ylim([-2,2])
colorbarpwn(0,20,'log','colorN',[1,1,1],'colorW',[0.2,0,0.8],'colorP',[1,0,1],'full', 5);
set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in');


%% 11. Force and flow snapshots aT4 = 0.33, h = 0.7c

clear;

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=33_p3=75_h3=0.7_ph=60_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.33_p3=75_h3=0.7c_ph=60.mat');

i = 0;

for n = [2,5,8,11,14]%,16]%1:3:16
    i = i + 1;
    x_60(:,:,i) = compiled_data(n).x;
    y_60(:,:,i) = compiled_data(n).y;
    foil2_coords_60(:,i) = compiled_data(n).foil2_coords;
    foil3_coords_60(:,i) = compiled_data(n).foil3_coords;
    u_60(:,:,i) = compiled_data(n).u;
    v_60(:,:,i) = compiled_data(n).v;
    vort = compiled_data(n).vort*0.061/0.33;
    isValid = compiled_data(n).isValid;
    vort(isValid == 0) = NaN; % foil mask
    vort_60(:,:,i) = vort;
    timeIndex_60(i) = compiled_data(n).timeStep + 1;
end

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=33_p3=75_h3=0.7_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.33_p3=75_h3=0.7c_ph=-120.mat');

i = 0;

for n = [14,17,20,23,3]%,4]%13:2:23
    i = i + 1;
    x_120(:,:,i) = compiled_data(n).x;
    y_120(:,:,i) = compiled_data(n).y;
    foil2_coords_120(:,i) = compiled_data(n).foil2_coords;
    foil3_coords_120(:,i) = compiled_data(n).foil3_coords;
    u_120(:,:,i) = compiled_data(n).u;
    v_120(:,:,i) = compiled_data(n).v;
    vort = compiled_data(n).vort*0.061/0.33;
    isValid = compiled_data(n).isValid;
    vort(isValid == 0) = NaN; % foil mask
    vort_120(:,:,i) = vort;
    timeIndex_120(i) = compiled_data(n).timeStep + 1;
end

timeIndex = mean(timeIndex_60,timeIndex_120);

load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221006_TandemPIV_3alphaRegimes\20221006_alpha=33_p3=75_h3=0.7_ph=60_A3E.mat');
[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, samplerate, transientcycs, foil_separation, flume_height);
out = filter_motor_noise_gromit(out, freq, samplerate, 20);
out = filter_motor_noise_wallace(out, freq, samplerate, 20);
res = calculate_forces(par, kin, out);
[     ~,          ~, hv_cyc3_60] = cycle_avg_data(kin.h3_comm/par.H3, kin.h3_vel, samplerate, 1); % h3 velocity
[     ~,          ~, CL3_cyc_60] = cycle_avg_data(kin.h3_comm/par.H3, res.CL3, samplerate, 1); % LiftC3
[toverT_60, heave_cyc3_60, CP3_cyc_60] = cycle_avg_data(kin.h3_comm/par.H3, (res.CPH3+res.CPP3), samplerate, 1); % PowerC3

load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221006_TandemPIV_3alphaRegimes\20221006_alpha=33_p3=75_h3=0.7_ph=-120_A3E.mat');
[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, samplerate, transientcycs, foil_separation, flume_height);
out = filter_motor_noise_gromit(out, freq, samplerate, 20);
out = filter_motor_noise_wallace(out, freq, samplerate, 20);
res = calculate_forces(par, kin, out);
[     ~,          ~, hv_cyc3_120] = cycle_avg_data(kin.h3_comm/par.H3, kin.h3_vel, samplerate, 1); % h3 velocity
[     ~,          ~, CL3_cyc_120] = cycle_avg_data(kin.h3_comm/par.H3, res.CL3, samplerate, 1); % LiftC3
[toverT_120, heave_cyc3_120, CP3_cyc_120] = cycle_avg_data(kin.h3_comm/par.H3, (res.CPH3+res.CPP3), samplerate, 1); % PowerC3



color_pitch = [0.4,0.4,0.4];
color_60 = [0.8,0,0];
color_120 = [0,0.5,0];
snapshot = {'a','b','c','d','e'};
yposition = [1.1, 1.1, 1.1, 1.1, 1.1];

figure(7)

for n = 1:5
    subplot(3,5,n)
    contourf(x_60(:,:,n),y_60(:,:,n),vort_60(:,:,n),'LineStyle','none','LevelStep',0.1);
    colormap(brewermap([],"-RdBu")); hold on;
    plot(foil2_coords_60([1,3,5],n), foil2_coords_60([2,4,6],n), 'k', 'LineWidth', 5); % plot leading foil
    plot(foil3_coords_60([1,3,5],n), foil3_coords_60([2,4,6],n), 'k', 'LineWidth', 5); hold off; % plot trailing foil
    axis equal
    caxis([-1,1]);
    xlim([0,4]);
    ylim([-2,2]);
    if n == 1
        set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','XTick',[],'XColor',color_60,'YColor',color_60);
        ylabel('$y/c$','interpreter','latex');
    else
        set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','XTick',[],'YTick',[],'XColor',color_60,'YColor',color_60);
    end
    subtitle(snapshot(n),'interpreter','latex','color',color_60);
end
clrbr = colorbar('Location','manual','units','normalized','Position',[0.93,0.71,0.01,0.21]);
clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
clrbr.Label.Interpreter = 'latex';

cyc_end = length(toverT_60);
cyc_range = 1:cyc_end;%*4/5;

subplot(3,5,[6,7]);
colororder({'b','k'})
yyaxis left
plot(toverT_60(cyc_range),mean(CL3_cyc_60(:,cyc_range)),'-','LineWidth',3,'color',color_60); hold on;
plot(toverT_60(timeIndex),mean(CL3_cyc_60(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_60);
plot(toverT_120(cyc_range),mean(CL3_cyc_120(:,cyc_range)),'-','LineWidth',3,'color',color_120);
plot(toverT_120(timeIndex),mean(CL3_cyc_120(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_120); hold off;
ylabel('$C_{\mathrm{L,tr}}$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
ylim([-5,5]);
yyaxis right
% plot(toverT_60(cyc_range),mean(heave_cyc3_60(:,cyc_range)),'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
% plot(toverT_60(timeIndex),mean(heave_cyc3_60(:,timeIndex)),'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
plot(toverT_60(cyc_range),mean(hv_cyc3_60(:,cyc_range))/U,'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
plot(toverT_60(timeIndex),mean(hv_cyc3_60(:,timeIndex))/U,'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
xline(toverT_60(timeIndex),'color','k','Alpha',0.2);
text(toverT_60(timeIndex), yposition, snapshot,'interpreter','latex','fontsize',17,'color',color_pitch); hold off;
ylim([-1,1]);
% ylabel('$h_{\mathrm{tr}}/c$','interpreter','latex');
ylabel('$\dot{h}_{\mathrm{tr}}/U_{\infty}$','interpreter','latex');
legend('$\psi_{1-2} = 60^{\circ}$,','','$\psi_{1-2} = -120^{\circ}$,','','$h^*_{\mathrm{tr}}$','interpreter','latex','location','southeast','orientation','horizontal','fontsize',18);
set(gca,'fontname','Times New Roman','fontsize',16,'linewidth',1.5,'TickDir','in');
box on;

subplot(3,5,[9,10]);
colororder({'b','k'})
yyaxis left
plot(toverT_60(cyc_range),mean(CP3_cyc_60(:,cyc_range)),'-','LineWidth',3,'color',color_60); hold on;
plot(toverT_60(timeIndex),mean(CP3_cyc_60(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_60);
plot(toverT_120(cyc_range),mean(CP3_cyc_120(:,cyc_range)),'-','LineWidth',3,'color',color_120);
plot(toverT_120(timeIndex),mean(CP3_cyc_120(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_120); hold off;
ylabel('$C_{\mathrm{P,tr}}$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
ylim([-1,1]);
yyaxis right
% plot(toverT_60(cyc_range),mean(heave_cyc3_60(:,cyc_range)),'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
% plot(toverT_60(timeIndex),mean(heave_cyc3_60(:,timeIndex)),'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
plot(toverT_60(cyc_range),mean(hv_cyc3_60(:,cyc_range))/U,'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
plot(toverT_60(timeIndex),mean(hv_cyc3_60(:,timeIndex))/U,'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
xline(toverT_60(timeIndex),'color','k','Alpha',0.2);
text(toverT_60(timeIndex), yposition, snapshot,'interpreter','latex','fontsize',17,'color',color_pitch); hold off;
ylim([-1,1]);
% ylabel('$h_{\mathrm{tr}}/c$','interpreter','latex');
ylabel('$\dot{h}_{\mathrm{tr}}/U_{\infty}$','interpreter','latex');
legend('$\psi_{1-2} = 60^{\circ}$,','','$\psi_{1-2} = -120^{\circ}$,','','$h^*_{\mathrm{tr}}$','interpreter','latex','location','southeast','orientation','horizontal','fontsize',18);
set(gca,'fontname','Times New Roman','fontsize',16,'linewidth',1.5,'TickDir','in');
box on;


for n = 1:5
    subplot(3,5,n+10)
    contourf(x_120(:,:,n),y_120(:,:,n),vort_120(:,:,n),'LineStyle','none','LevelStep',0.1);
    colormap(brewermap([],"-RdBu")); hold on;
    plot(foil2_coords_120([1,3,5],n), foil2_coords_120([2,4,6],n), 'k', 'LineWidth', 5); % plot leading foil
    plot(foil3_coords_120([1,3,5],n), foil3_coords_120([2,4,6],n), 'k', 'LineWidth', 5); hold off; % plot trailing foil
    axis equal
    caxis([-1,1]);
    xlim([0,4]);
    ylim([-2,2]);
    if n == 1
        set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','XColor',color_120,'YColor',color_120);
        ylabel('$y/c$','interpreter','latex');
        xlabel('$x/c$','interpreter','latex');
    else
        set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','YTick',[],'XColor',color_120,'YColor',color_120);
        xlabel('$x/c$','interpreter','latex');
    end
    subtitle(snapshot(n),'interpreter','latex','color',color_120);
end
clrbr = colorbar('Location','manual','units','normalized','Position',[0.93,0.11,0.01,0.21]);
clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
clrbr.Label.Interpreter = 'latex';

disp('done');

%% 12. Force and flow snapshots aT4 = 0.68, h = 0.7c

clear;

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=68_p3=75_h3=0.7_ph=60_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.68_p3=75_h3=0.7c_ph=60.mat');

i = 0;

for n = [2,5,8,11,14]%,16]%1:3:16
    i = i + 1;
    x_60(:,:,i) = compiled_data(n).x;
    y_60(:,:,i) = compiled_data(n).y;
    foil2_coords_60(:,i) = compiled_data(n).foil2_coords;
    foil3_coords_60(:,i) = compiled_data(n).foil3_coords;
    u_60(:,:,i) = compiled_data(n).u;
    v_60(:,:,i) = compiled_data(n).v;
    vort = compiled_data(n).vort*0.061/0.33;
    isValid = compiled_data(n).isValid;
    vort(isValid == 0) = NaN; % foil mask
    vort_60(:,:,i) = vort;
    timeIndex_60(i) = compiled_data(n).timeStep + 1;
end

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=68_p3=75_h3=0.7_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.68_p3=75_h3=0.7c_ph=-120.mat');

i = 0;

for n = [2,5,8,11,14]
    i = i + 1;
    x_120(:,:,i) = compiled_data(n).x;
    y_120(:,:,i) = compiled_data(n).y;
    foil2_coords_120(:,i) = compiled_data(n).foil2_coords;
    foil3_coords_120(:,i) = compiled_data(n).foil3_coords;
    u_120(:,:,i) = compiled_data(n).u;
    v_120(:,:,i) = compiled_data(n).v;
    vort = compiled_data(n).vort*0.061/0.33;
    isValid = compiled_data(n).isValid;
    vort(isValid == 0) = NaN; % foil mask
    vort_120(:,:,i) = vort;
    timeIndex_120(i) = compiled_data(n).timeStep + 1;
end

timeIndex = mean(timeIndex_60,timeIndex_120);

load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221006_TandemPIV_3alphaRegimes\20221006_alpha=68_p3=75_h3=0.7_ph=60_A3E.mat');
[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, samplerate, transientcycs, foil_separation, flume_height);
out = filter_motor_noise_gromit(out, freq, samplerate, 20);
out = filter_motor_noise_wallace(out, freq, samplerate, 20);
res = calculate_forces(par, kin, out);
[     ~,          ~, hv_cyc3_60] = cycle_avg_data(kin.h3_comm/par.H3, kin.h3_vel, samplerate, 1); % h3 velocity
[     ~,          ~, CL3_cyc_60] = cycle_avg_data(kin.h3_comm/par.H3, res.CL3, samplerate, 1); % LiftC3
[toverT_60, heave_cyc3_60, CP3_cyc_60] = cycle_avg_data(kin.h3_comm/par.H3, (res.CPH3+res.CPP3), samplerate, 1); % PowerC3

load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221006_TandemPIV_3alphaRegimes\20221006_alpha=68_p3=75_h3=0.7_ph=-120_A3E.mat');
[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, samplerate, transientcycs, foil_separation, flume_height);
out = filter_motor_noise_gromit(out, freq, samplerate, 20);
out = filter_motor_noise_wallace(out, freq, samplerate, 20);
res = calculate_forces(par, kin, out);
[     ~,          ~, hv_cyc3_120] = cycle_avg_data(kin.h3_comm/par.H3, kin.h3_vel, samplerate, 1); % h3 velocity
[     ~,          ~, CL3_cyc_120] = cycle_avg_data(kin.h3_comm/par.H3, res.CL3, samplerate, 1); % LiftC3
[toverT_120, heave_cyc3_120, CP3_cyc_120] = cycle_avg_data(kin.h3_comm/par.H3, (res.CPH3+res.CPP3), samplerate, 1); % PowerC3


color_pitch = [0.4,0.4,0.4];
color_60 = [0.8,0,0];
color_120 = [0,0.5,0];
snapshot = {'a','b','c','d','e'};
yposition = [1.1, 1.1, 1.1, 1.1, 1.1];

figure(7)

for n = 1:5
    subplot(3,5,n)
    contourf(x_60(:,:,n),y_60(:,:,n),vort_60(:,:,n),'LineStyle','none','LevelStep',0.1);
    colormap(brewermap([],"-RdBu")); hold on;
    plot(foil2_coords_60([1,3,5],n), foil2_coords_60([2,4,6],n), 'k', 'LineWidth', 5); % plot leading foil
    plot(foil3_coords_60([1,3,5],n), foil3_coords_60([2,4,6],n), 'k', 'LineWidth', 5); hold off; % plot trailing foil
    axis equal
    caxis([-1,1]);
    xlim([0,4]);
    ylim([-2,2]);
    if n == 1
        set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','XTick',[],'XColor',color_60,'YColor',color_60);
        ylabel('$y/c$','interpreter','latex');
    else
        set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','XTick',[],'YTick',[],'XColor',color_60,'YColor',color_60);
    end
    subtitle(snapshot(n),'interpreter','latex','color',color_60);
end
clrbr = colorbar('Location','manual','units','normalized','Position',[0.93,0.71,0.01,0.21]);
clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
clrbr.Label.Interpreter = 'latex';

cyc_end = length(toverT_60);
cyc_range = 1:cyc_end;%*4/5;

subplot(3,5,[6,7]);
colororder({'b','k'})
yyaxis left
plot(toverT_60(cyc_range),mean(CL3_cyc_60(:,cyc_range)),'-','LineWidth',3,'color',color_60); hold on;
plot(toverT_60(timeIndex),mean(CL3_cyc_60(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_60);
plot(toverT_120(cyc_range),mean(CL3_cyc_120(:,cyc_range)),'-','LineWidth',3,'color',color_120);
plot(toverT_120(timeIndex),mean(CL3_cyc_120(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_120); hold off;
ylabel('$C_{\mathrm{L,tr}}$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
ylim([-5,5]);
yyaxis right
% plot(toverT_60(cyc_range),mean(heave_cyc3_60(:,cyc_range)),'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
% plot(toverT_60(timeIndex),mean(heave_cyc3_60(:,timeIndex)),'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
plot(toverT_60(cyc_range),mean(hv_cyc3_60(:,cyc_range))/U,'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
plot(toverT_60(timeIndex),mean(hv_cyc3_60(:,timeIndex))/U,'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
xline(toverT_60(timeIndex),'color','k','Alpha',0.2);
text(toverT_60(timeIndex), yposition, snapshot,'interpreter','latex','fontsize',17,'color',color_pitch); hold off;
ylim([-1,1]);
% ylabel('$h_{\mathrm{tr}}/c$','interpreter','latex');
ylabel('$\dot{h}_{\mathrm{tr}}/U_{\infty}$','interpreter','latex');
legend('$\psi_{1-2} = 60^{\circ}$,','','$\psi_{1-2} = -120^{\circ}$,','','$h^*_{\mathrm{tr}}$','interpreter','latex','location','southeast','orientation','horizontal','fontsize',18);
set(gca,'fontname','Times New Roman','fontsize',16,'linewidth',1.5,'TickDir','in');
box on;

subplot(3,5,[9,10]);
colororder({'b','k'})
yyaxis left
plot(toverT_60(cyc_range),mean(CP3_cyc_60(:,cyc_range)),'-','LineWidth',3,'color',color_60); hold on;
plot(toverT_60(timeIndex),mean(CP3_cyc_60(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_60);
plot(toverT_120(cyc_range),mean(CP3_cyc_120(:,cyc_range)),'-','LineWidth',3,'color',color_120);
plot(toverT_120(timeIndex),mean(CP3_cyc_120(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_120); hold off;
ylabel('$C_{\mathrm{P,tr}}$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
ylim([-1,1]);
yyaxis right
% plot(toverT_60(cyc_range),mean(heave_cyc3_60(:,cyc_range)),'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
% plot(toverT_60(timeIndex),mean(heave_cyc3_60(:,timeIndex)),'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
plot(toverT_60(cyc_range),mean(hv_cyc3_60(:,cyc_range))/U,'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
plot(toverT_60(timeIndex),mean(hv_cyc3_60(:,timeIndex))/U,'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
xline(toverT_60(timeIndex),'color','k','Alpha',0.2);
text(toverT_60(timeIndex), yposition, snapshot,'interpreter','latex','fontsize',17,'color',color_pitch); hold off;
ylim([-1,1]);
% ylabel('$h_{\mathrm{tr}}/c$','interpreter','latex');
ylabel('$\dot{h}_{\mathrm{tr}}/U_{\infty}$','interpreter','latex');
legend('$\psi_{1-2} = 60^{\circ}$,','','$\psi_{1-2} = -120^{\circ}$,','','$h^*_{\mathrm{tr}}$','interpreter','latex','location','southeast','orientation','horizontal','fontsize',18);
set(gca,'fontname','Times New Roman','fontsize',16,'linewidth',1.5,'TickDir','in');
box on;


for n = 1:5
    subplot(3,5,n+10)
    contourf(x_120(:,:,n),y_120(:,:,n),vort_120(:,:,n),'LineStyle','none','LevelStep',0.1);
    colormap(brewermap([],"-RdBu")); hold on;
    plot(foil2_coords_120([1,3,5],n), foil2_coords_120([2,4,6],n), 'k', 'LineWidth', 5); % plot leading foil
    plot(foil3_coords_120([1,3,5],n), foil3_coords_120([2,4,6],n), 'k', 'LineWidth', 5); hold off; % plot trailing foil
    axis equal
    caxis([-1,1]);
    xlim([0,4]);
    ylim([-2,2]);
    if n == 1
        set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','XColor',color_120,'YColor',color_120);
        ylabel('$y/c$','interpreter','latex');
        xlabel('$x/c$','interpreter','latex');
    else
        set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','YTick',[],'XColor',color_120,'YColor',color_120);
        xlabel('$x/c$','interpreter','latex');
    end
    subtitle(snapshot(n),'interpreter','latex','color',color_120);
end
clrbr = colorbar('Location','manual','units','normalized','Position',[0.93,0.11,0.01,0.21]);
clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
clrbr.Label.Interpreter = 'latex';

disp('done');


%% 12. Force and flow snapshots aT4 = 0.68, h = 0.8c

clear;

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=68_p3=75_h3=0.8_ph=60_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.68_p3=75_h3=0.8c_ph=60.mat');

i = 0;

for n = [2,5,8,11,14]%,16]%1:3:16
    i = i + 1;
    x_60(:,:,i) = compiled_data(n).x;
    y_60(:,:,i) = compiled_data(n).y;
    foil2_coords_60(:,i) = compiled_data(n).foil2_coords;
    foil3_coords_60(:,i) = compiled_data(n).foil3_coords;
    u_60(:,:,i) = compiled_data(n).u;
    v_60(:,:,i) = compiled_data(n).v;
    vort = compiled_data(n).vort*0.061/0.33;
    isValid = compiled_data(n).isValid;
    vort(isValid == 0) = NaN; % foil mask
    vort_60(:,:,i) = vort;
    timeIndex_60(i) = compiled_data(n).timeStep + 1;
end

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=68_p3=75_h3=0.8_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.68_p3=75_h3=0.8c_ph=-120.mat');

i = 0;

for n = [14,17,20,23,3]
    i = i + 1;
    x_120(:,:,i) = compiled_data(n).x;
    y_120(:,:,i) = compiled_data(n).y;
    foil2_coords_120(:,i) = compiled_data(n).foil2_coords;
    foil3_coords_120(:,i) = compiled_data(n).foil3_coords;
    u_120(:,:,i) = compiled_data(n).u;
    v_120(:,:,i) = compiled_data(n).v;
    vort = compiled_data(n).vort*0.061/0.33;
    isValid = compiled_data(n).isValid;
    vort(isValid == 0) = NaN; % foil mask
    vort_120(:,:,i) = vort;
    timeIndex_120(i) = compiled_data(n).timeStep + 1;
end

timeIndex = mean(timeIndex_60,timeIndex_120);

load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221006_TandemPIV_3alphaRegimes\20221006_alpha=68_p3=75_h3=0.8_ph=60_A3E.mat');
[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, samplerate, transientcycs, foil_separation, flume_height);
out = filter_motor_noise_gromit(out, freq, samplerate, 20);
out = filter_motor_noise_wallace(out, freq, samplerate, 20);
res = calculate_forces(par, kin, out);
[     ~,          ~, hv_cyc3_60] = cycle_avg_data(kin.h3_comm/par.H3, kin.h3_vel, samplerate, 1); % h3 velocity
[     ~,          ~, CL3_cyc_60] = cycle_avg_data(kin.h3_comm/par.H3, res.CL3, samplerate, 1); % LiftC3
[toverT_60, heave_cyc3_60, CP3_cyc_60] = cycle_avg_data(kin.h3_comm/par.H3, (res.CPH3+res.CPP3), samplerate, 1); % PowerC3

load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221006_TandemPIV_3alphaRegimes\20221006_alpha=68_p3=75_h3=0.8_ph=-120_A3E.mat');
[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, samplerate, transientcycs, foil_separation, flume_height);
out = filter_motor_noise_gromit(out, freq, samplerate, 20);
out = filter_motor_noise_wallace(out, freq, samplerate, 20);
res = calculate_forces(par, kin, out);
[     ~,          ~, hv_cyc3_120] = cycle_avg_data(kin.h3_comm/par.H3, kin.h3_vel, samplerate, 1); % h3 velocity
[     ~,          ~, CL3_cyc_120] = cycle_avg_data(kin.h3_comm/par.H3, res.CL3, samplerate, 1); % LiftC3
[toverT_120, heave_cyc3_120, CP3_cyc_120] = cycle_avg_data(kin.h3_comm/par.H3, (res.CPH3+res.CPP3), samplerate, 1); % PowerC3


color_pitch = [0.4,0.4,0.4];
color_60 = [0.8,0,0];
color_120 = [0,0.5,0];
snapshot = {'a','b','c','d','e'};
yposition = [1.1, 1.1, 1.1, 1.1, 1.1];

figure(7)

for n = 1:5
    subplot(3,5,n)
    contourf(x_60(:,:,n),y_60(:,:,n),vort_60(:,:,n),'LineStyle','none','LevelStep',0.1);
    colormap(brewermap([],"-RdBu")); hold on;
    plot(foil2_coords_60([1,3,5],n), foil2_coords_60([2,4,6],n), 'k', 'LineWidth', 5); % plot leading foil
    plot(foil3_coords_60([1,3,5],n), foil3_coords_60([2,4,6],n), 'k', 'LineWidth', 5); hold off; % plot trailing foil
    axis equal
    caxis([-1,1]);
    xlim([0,4]);
    ylim([-2,2]);
    if n == 1
        set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','XTick',[],'XColor',color_60,'YColor',color_60);
        ylabel('$y/c$','interpreter','latex');
    else
        set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','XTick',[],'YTick',[],'XColor',color_60,'YColor',color_60);
    end
    subtitle(snapshot(n),'interpreter','latex','color',color_60);
end
clrbr = colorbar('Location','manual','units','normalized','Position',[0.93,0.71,0.01,0.21]);
clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
clrbr.Label.Interpreter = 'latex';

cyc_end = length(toverT_60);
cyc_range = 1:cyc_end;%*4/5;

subplot(3,5,[6,7]);
colororder({'b','k'})
yyaxis left
plot(toverT_60(cyc_range),mean(CL3_cyc_60(:,cyc_range)),'-','LineWidth',3,'color',color_60); hold on;
plot(toverT_60(timeIndex),mean(CL3_cyc_60(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_60);
plot(toverT_120(cyc_range),mean(CL3_cyc_120(:,cyc_range)),'-','LineWidth',3,'color',color_120);
plot(toverT_120(timeIndex),mean(CL3_cyc_120(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_120); hold off;
ylabel('$C_{\mathrm{L,tr}}$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
ylim([-5,5]);
yyaxis right
% plot(toverT_60(cyc_range),mean(heave_cyc3_60(:,cyc_range)),'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
% plot(toverT_60(timeIndex),mean(heave_cyc3_60(:,timeIndex)),'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
plot(toverT_60(cyc_range),mean(hv_cyc3_60(:,cyc_range))/U,'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
plot(toverT_60(timeIndex),mean(hv_cyc3_60(:,timeIndex))/U,'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
xline(toverT_60(timeIndex),'color','k','Alpha',0.2);
text(toverT_60(timeIndex), yposition, snapshot,'interpreter','latex','fontsize',17,'color',color_pitch); hold off;
ylim([-1,1]);
% ylabel('$h_{\mathrm{tr}}/c$','interpreter','latex');
ylabel('$\dot{h}_{\mathrm{tr}}/U_{\infty}$','interpreter','latex');
legend('$\psi_{1-2} = 60^{\circ}$,','','$\psi_{1-2} = -120^{\circ}$,','','$h^*_{\mathrm{tr}}$','interpreter','latex','location','southeast','orientation','horizontal','fontsize',18);
set(gca,'fontname','Times New Roman','fontsize',16,'linewidth',1.5,'TickDir','in');
box on;

subplot(3,5,[9,10]);
colororder({'b','k'})
yyaxis left
plot(toverT_60(cyc_range),mean(CP3_cyc_60(:,cyc_range)),'-','LineWidth',3,'color',color_60); hold on;
plot(toverT_60(timeIndex),mean(CP3_cyc_60(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_60);
plot(toverT_120(cyc_range),mean(CP3_cyc_120(:,cyc_range)),'-','LineWidth',3,'color',color_120);
plot(toverT_120(timeIndex),mean(CP3_cyc_120(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_120); hold off;
ylabel('$C_{\mathrm{P,tr}}$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
ylim([-1,1]);
yyaxis right
% plot(toverT_60(cyc_range),mean(heave_cyc3_60(:,cyc_range)),'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
% plot(toverT_60(timeIndex),mean(heave_cyc3_60(:,timeIndex)),'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
plot(toverT_60(cyc_range),mean(hv_cyc3_60(:,cyc_range))/U,'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
plot(toverT_60(timeIndex),mean(hv_cyc3_60(:,timeIndex))/U,'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
xline(toverT_60(timeIndex),'color','k','Alpha',0.2);
text(toverT_60(timeIndex), yposition, snapshot,'interpreter','latex','fontsize',17,'color',color_pitch); hold off;
ylim([-1,1]);
% ylabel('$h_{\mathrm{tr}}/c$','interpreter','latex');
ylabel('$\dot{h}_{\mathrm{tr}}/U_{\infty}$','interpreter','latex');
legend('$\psi_{1-2} = 60^{\circ}$,','','$\psi_{1-2} = -120^{\circ}$,','','$h^*_{\mathrm{tr}}$','interpreter','latex','location','southeast','orientation','horizontal','fontsize',18);
set(gca,'fontname','Times New Roman','fontsize',16,'linewidth',1.5,'TickDir','in');
box on;


for n = 1:5
    subplot(3,5,n+10)
    contourf(x_120(:,:,n),y_120(:,:,n),vort_120(:,:,n),'LineStyle','none','LevelStep',0.1);
    colormap(brewermap([],"-RdBu")); hold on;
    plot(foil2_coords_120([1,3,5],n), foil2_coords_120([2,4,6],n), 'k', 'LineWidth', 5); % plot leading foil
    plot(foil3_coords_120([1,3,5],n), foil3_coords_120([2,4,6],n), 'k', 'LineWidth', 5); hold off; % plot trailing foil
    axis equal
    caxis([-1,1]);
    xlim([0,4]);
    ylim([-2,2]);
    if n == 1
        set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','XColor',color_120,'YColor',color_120);
        ylabel('$y/c$','interpreter','latex');
        xlabel('$x/c$','interpreter','latex');
    else
        set(gca,'fontname','Times New Roman','fontsize',18,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','YTick',[],'XColor',color_120,'YColor',color_120);
        xlabel('$x/c$','interpreter','latex');
    end
    subtitle(snapshot(n),'interpreter','latex','color',color_120);
end
clrbr = colorbar('Location','manual','units','normalized','Position',[0.93,0.11,0.01,0.21]);
clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
clrbr.Label.Interpreter = 'latex';

disp('done');


%% Mean wake velocity at a = 0.68


clear;

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=68_p3=75_h3=0.7_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.68_p3=75_h3=0.7c_ph=-120.mat');

frame = 23;

for n = 1:frame
    x(:,:,n) = compiled_data(n).x;
    y(:,:,n) = compiled_data(n).y;
    foil2_coords(:,n) = compiled_data(n).foil2_coords;
    foil3_coords(:,n) = compiled_data(n).foil3_coords;
    u(:,:,n) = compiled_data(n).u*0.33;
    v(:,:,n) = compiled_data(n).v*0.33;
    U_sm = sqrt((compiled_data(n).u).^2 + (compiled_data(n).v).^2);
    isValid = compiled_data(n).isValid;
    U_sm(isValid == 0) = NaN; % foil mask
    U(:,:,n) = U_sm;
    timeIndex(n) = compiled_data(n).timeStep + 1;
end

k = 2;

contourf(x(:,:,k),y(:,:,k),u(:,:,k)); hold on;
plot(x(226,123,k), y(226,123,k), 'ro', 'linewidth', 3, 'markersize', 3);
plot(x(226,123,k), y(226,123,k), 'ro', 'linewidth', 3, 'markersize', 15); hold off;

u_avg = mean(squeeze(U(226,123,:)));

%% 3. aT4 vs eta - Ribeiro et al data

clear;

% Load .csv file with data from Webplotdigitizer
filename = '\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221028_RibeiroEtAl_2022_aT4_vs_eff\Default Dataset.csv';
A = readmatrix(filename);


figure(5)

plot(deg2rad(A(:,1)),A(:,2),'ko','markersize',9,'linewidth',1.5); hold on;
xline(deg2rad([11.7,29.3]),'b','linewidth',2); hold off;
ylim([0,0.4]);
xlim([0,1]);
xlabel('$\alpha_{T/4}[rad]$','interpreter','latex');
ylabel('$\eta$','interpreter','latex','rotation',0);
xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]);
yticks([0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40]);
set(gca, 'fontname','Times New Roman','fontsize',20,'LineWidth',1.5,'LabelFontSizeMultiplier',1.5);
grid on;
box on;
