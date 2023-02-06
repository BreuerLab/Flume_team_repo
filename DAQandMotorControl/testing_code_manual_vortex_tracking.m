%% testing code - manual vortex tracking
% 20221213

clear;

frame = 23;
half = 11;

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=33_p3=75_h3=0.7_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.33_p3=75_h3=0.7c_ph=-120.mat');

i = 0;

for n = [12:23,1:9]
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

for n = [1:21]%[12:23,1:9]
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
    for m = 1:21

        figure('name', 'a33', 'WindowState', 'maximized');
        contourf(x_33(:,:,m), y_33(:,:,m), Q_33(:,:,m),'LineStyle','none');
        axis equal
        caxis([0,20])
        xlim([-3.5,3.5])
        ylim([-2,2])
    %     colorbarpwn(0,20,'log','colorN',[1,1,1],'colorW',[0.2,0,0.8],'colorP',[1,0,1],'full', 5);
        LEV_ctr_33 = drawpoint();
        LEV_coords_33(do_it_again,m,:) = LEV_ctr_33.Position;
        close 'a33';

        figure('name', 'a68', 'WindowState', 'maximized');
        contourf(x_68(:,:,m), y_68(:,:,m), Q_68(:,:,m),'LineStyle','none');
        axis equal
        caxis([0,20])
        xlim([-3.5,3.5])
        ylim([-2,2])
    %     colorbarpwn(0,20,'log','colorN',[1,1,1],'colorW',[0.2,0,0.8],'colorP',[1,0,1],'full', 5);
        LEV_ctr_68 = drawpoint();
        LEV_coords_68(do_it_again,m,:) = LEV_ctr_68.Position;
        close 'a68';

    end
end

LEV_coords_33_avg = squeeze(mean(LEV_coords_33,1));
LEV_coords_68_avg = squeeze(mean(LEV_coords_68,1));

save('20231217_Main_vortex_tracking_a=0.33_p3=75_h3=0.7c_ph=-120.mat','LEV_coords_33','LEV_coords_33_avg');
save('20231217_Main_vortex_tracking_a=0.68_p3=75_h3=0.7c_ph=-120.mat','LEV_coords_68','LEV_coords_68_avg');
%%
load('20231217_Main_vortex_tracking_a=0.33_p3=75_h3=0.7c_ph=-120.mat');
load('20231217_Main_vortex_tracking_a=0.68_p3=75_h3=0.7c_ph=-120.mat');

LEV_coords_33_avg(:,1) = LEV_coords_33_avg(:,1) + 3;
LEV_coords_68_avg(:,1) = LEV_coords_68_avg(:,1) + 3;

figure(1)

toime = 1:half;
plot(LEV_coords_33_avg(:,1),LEV_coords_33_avg(:,2),'o','color',[0.4940 0.1840 0.5560],'MarkerSize',10,'LineWidth',2,'MarkerFaceColor',[0.4940 0.1840 0.5560]); hold on;
plot(LEV_coords_68_avg(:,1),LEV_coords_68_avg(:,2),'s','color',[0.3010 0.7450 0.9330],'MarkerSize',10,'LineWidth',2,'MarkerFaceColor',[0.3010 0.7450 0.9330]); hold off;
legend('$\alpha_{T/4} = 0.33$','$\alpha_{T/4} = 0.68$','interpreter','latex','fontsize',16,'location','southeast');
% tit = 'Primary Vortex Trajectory';
xlim([0,6])
% title(tit,'interpreter','latex','fontsize',24);
set(gca,'fontname','Times New Roman','fontsize',16,'linewidth',1.5,'TickDir','in');
xlabel('$x/c$','interpreter','latex');
ylabel('$y/c$','interpreter','latex');
grid on;
% axis equal
