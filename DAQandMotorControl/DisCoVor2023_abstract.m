%% DisCorVor 2023 - Abstract Figures

%% 13. Force and flow snapshots aT4 = 0.33, h = 0.7c

clear;

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=33_p3=75_h3=0.7_ph=60_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export\TandemFoil_aT4=0.33_p3=75_h3=0.7c_ph=60.mat');

i = 0;

for n = [5,11,14]%,16]%1:3:16
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

for n = [17,23,3]
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

load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221114_TandemMonday_redux_LeadingAlphaSweep\data\20221114_TandemMonday_leadingAlphaSweep_PHPh_aT4=0.33rad,p3=75deg,h3=0.7c,ph=60.mat');
[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, samplerate, transientcycs, foil_separation, flume_height);
out = filter_motor_noise_gromit(out, freq, samplerate, 15);
out = filter_motor_noise_wallace(out, freq, samplerate, 15);
res = calculate_forces(par, kin, out);
[     ~,          ~, hv_cyc3_60] = cycle_avg_data(kin.h3_comm/par.H3, kin.h3_vel, samplerate, 1); % h3 velocity
[     ~,          ~, CL3_cyc_60] = cycle_avg_data(kin.h3_comm/par.H3, res.CL3, samplerate, 1); % LiftC3
[toverT_60, heave_cyc3_60, CP3_cyc_60] = cycle_avg_data(kin.h3_comm/par.H3, (res.CPH3+res.CPP3), samplerate, 1); % PowerC3

load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221114_TandemMonday_redux_LeadingAlphaSweep\data\20221114_TandemMonday_leadingAlphaSweep_PHPh_aT4=0.33rad,p3=75deg,h3=0.7c,ph=-120.mat');
[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, samplerate, transientcycs, foil_separation, flume_height);
out = filter_motor_noise_gromit(out, freq, samplerate, 15);
out = filter_motor_noise_wallace(out, freq, samplerate, 15);
res = calculate_forces(par, kin, out);
[     ~,          ~, hv_cyc3_120] = cycle_avg_data(kin.h3_comm/par.H3, kin.h3_vel, samplerate, 1); % h3 velocity
[     ~,          ~, CL3_cyc_120] = cycle_avg_data(kin.h3_comm/par.H3, res.CL3, samplerate, 1); % LiftC3
[toverT_120, heave_cyc3_120, CP3_cyc_120] = cycle_avg_data(kin.h3_comm/par.H3, (res.CPH3+res.CPP3), samplerate, 1); % PowerC3

%%
color_pitch = [0.4,0.4,0.4];
color_60 = [0,0.4470,0.7410];
color_120 = [0.6350,0.0780,0.1840];
snapshot = {'\textbf{a}','\textbf{b}','\textbf{c}'};
yposition = [0.85, 0.85, 0.85];
fz = 16; % font size

figure(1)

ii = 1:2:5;
for n = 1:3
    subplot(8,6,[ii(n),ii(n)+1,ii(n)+6,ii(n)+7])
    contourf(x_60(:,:,n),y_60(:,:,n),vort_60(:,:,n),'LineStyle','none','LevelStep',0.1);
    colormap(brewermap([],"-RdBu")); hold on;
    plot(foil2_coords_60([1,3,5],n), foil2_coords_60([2,4,6],n), 'k', 'LineWidth', 5); % plot leading foil
    plot(foil3_coords_60([1,3,5],n), foil3_coords_60([2,4,6],n), 'k', 'LineWidth', 5); hold off; % plot trailing foil
    axis equal
    caxis([-1,1]);
    xlim([0,4]);
    ylim([-2,2]);
    if n == 1
        set(gca,'fontname','Times New Roman','fontsize',fz,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','xticklabels',{},'XColor',color_60,'YColor',color_60);
        ylabel('$y/c$','interpreter','latex');
%         xlabel('$x/c$','interpreter','latex');
    else
        set(gca,'fontname','Times New Roman','fontsize',fz,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','xticklabels',{},'yticklabels',{},'XColor',color_60,'YColor',color_60);
%         xlabel('$x/c$','interpreter','latex');
    end
%     subtitle(snapshot(n),'interpreter','latex','color',color_60,'fontsize',24);
    text(3.2, -1.5, snapshot(n),'interpreter','latex','fontsize',fz+2,'color',color_60);
end
clrbr = colorbar('Location','manual','units','normalized','Position',[0.13,0.95,0.78,0.02],'LineWidth',1.5,'orientation','horizontal');
clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
clrbr.Label.Interpreter = 'latex';

cyc_end = length(toverT_60);
cyc_range = 1:cyc_end;%*4/5;

subplot(8,6,[13:1:24]);
colororder({'b','k'})
yyaxis left
plot(toverT_60(cyc_range),mean(CL3_cyc_60(:,cyc_range)),'-','LineWidth',3,'color',color_60); hold on;
plot(toverT_60(timeIndex),mean(CL3_cyc_60(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_60);
plot(toverT_120(cyc_range),mean(CL3_cyc_120(:,cyc_range)),'-','LineWidth',3,'color',color_120);
plot(toverT_120(timeIndex),mean(CL3_cyc_120(:,timeIndex)),'s','MarkerSize',12,'LineWidth',1.5,'color',color_120); hold off;
ylabel('$C_{\mathrm{L,tr}}$','interpreter','latex');
% xlabel('$t/T$','interpreter','latex');
ylim([-5,5]);
yyaxis right
% plot(toverT_60(cyc_range),mean(heave_cyc3_60(:,cyc_range)),'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
% plot(toverT_60(timeIndex),mean(heave_cyc3_60(:,timeIndex)),'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
plot(toverT_60(cyc_range),mean(hv_cyc3_60(:,cyc_range))/U,'LineWidth',1.5,'color',[0.5,0.5,0.5]); hold on;
plot(toverT_60(timeIndex),mean(hv_cyc3_60(:,timeIndex))/U,'s','MarkerSize',10,'LineWidth',1.5,'color',color_pitch);
xline(toverT_60(timeIndex),'color','k','Alpha',0.2);
text(toverT_60(timeIndex), yposition, snapshot,'interpreter','latex','fontsize',fz,'color',color_pitch); hold off;
ylim([-1,1]);
% ylabel('$h_{\mathrm{tr}}/c$','interpreter','latex');
ylabel('$\dot{h}_{\mathrm{tr}}/U_{\infty}$','interpreter','latex');
legend('$\psi_{1-2} = 60^{\circ}$,','','$\psi_{1-2} = -120^{\circ}$,','','$\dot{h}_{\mathrm{tr}}/U_{\infty}$','interpreter','latex','location','southwest','orientation','horizontal','fontsize',18);
set(gca,'fontname','Times New Roman','fontsize',fz,'linewidth',1.5,'TickDir','in','LabelFontSizeMultiplier',1.2);
box on;

subplot(8,6,[25:1:36]);
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
text(toverT_60(timeIndex), -yposition, snapshot,'interpreter','latex','fontsize',fz,'color',color_pitch); hold off;
ylim([-1,1]);
% ylabel('$h_{\mathrm{tr}}/c$','interpreter','latex');
ylabel('$\dot{h}_{\mathrm{tr}}/U_{\infty}$','interpreter','latex');
% legend('$\psi_{1-2} = 60^{\circ}$,','','$\psi_{1-2} = -120^{\circ}$,','','$h^*_{\mathrm{tr}}$','interpreter','latex','location','south','orientation','horizontal','fontsize',18);
set(gca,'fontname','Times New Roman','fontsize',fz,'linewidth',1.5,'TickDir','in','xticklabels',{},'LabelFontSizeMultiplier',1.2);
box on;


ii = 37:2:41;
for n = 1:3
    subplot(8,6,[ii(n),ii(n)+1,ii(n)+6,ii(n)+7])
    contourf(x_120(:,:,n),y_120(:,:,n),vort_120(:,:,n),'LineStyle','none','LevelStep',0.1);
    colormap(brewermap([],"-RdBu")); hold on;
    plot(foil2_coords_120([1,3,5],n), foil2_coords_120([2,4,6],n), 'k', 'LineWidth', 5); % plot leading foil
    plot(foil3_coords_120([1,3,5],n), foil3_coords_120([2,4,6],n), 'k', 'LineWidth', 5); hold off; % plot trailing foil
    axis equal
    caxis([-1,1]);
    xlim([0,4]);
    ylim([-2,2]);
    if n == 1
        set(gca,'fontname','Times New Roman','fontsize',fz,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','XColor',color_120,'YColor',color_120);
        ylabel('$y/c$','interpreter','latex');
        xlabel('$x/c$','interpreter','latex');
    else
        set(gca,'fontname','Times New Roman','fontsize',fz,'linewidth',1.5,'color',[0.8,0.8,0.8],'TickDir','in','YTick',[],'XColor',color_120,'YColor',color_120);
        xlabel('$x/c$','interpreter','latex');
    end
%     subtitle(snapshot(n),'interpreter','latex','color',color_120,'fontsize',24);
    text(3.2, -1.5, snapshot(n),'interpreter','latex','fontsize',fz+2,'color',color_120);
end
% clrbr = colorbar('Location','manual','units','normalized','Position',[0.93,0.155,0.01,0.25],'LineWidth',1.5);
% clrbr.Label.String = '$\omega_{z}c/U_{\infty}$';
% clrbr.Label.Interpreter = 'latex';

disp('done');