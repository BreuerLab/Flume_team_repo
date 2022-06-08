
%         % Plot heave spectrum
%     figure
%     hold on; grid on;
%     plot([f_heave_dom f_heave_dom],[-60 30],'Color','red');
%     plot(f_heave,10*log10(heave_powerspec));
% %     plot([f_force_dom f_force_dom],[-60 30],'Color','red','DisplayName','Frequency dominant');
% %     plot(f_force,10*log10(force_powerspec),'Color','black','DisplayName','Force spectrum'); % plot force spectrum
%     xlabel('Frequency (Hz)')
%     ylabel('Heave PSD (dB/Hz)')
%     xlim([0 10]) 
%     ylim([0 30])
% %     title(['CPS ring down, simulated: m = ',num2str(M),' kg, f_{natural} = ',num2str(f_nat),' Hz, \zeta = ',num2str(zeta)]);
%     hold off

    figure
    hold on; grid on;

    % Plot heave spectrum

%     plot([f_heave_dom f_heave_dom],[-60 30],'Color','red','DisplayName','Frequency dominant');
%     plot(f_heave,10*log10(heave_powerspec));
%     plot([f_force_dom f_force_dom],[-60 30],'Color','red');
    plot([1 1],[-60 60],'Color','red','LineStyle','--');
    plot([f_vortex/freq f_vortex/freq],[-60 60],'color','blue','LineStyle','--')
    plot([f_flume/freq,f_flume/freq],[-60,60],'color','cyan','LineStyle','--')
    plot(f_force/freq,10*log10(force_powerspec),'Color','black'); % plot force spectrum
    scatter(forcespec_peaklocs,forcespec_peakpowers,'Marker','o');
    xlabel('Frequency (f/f_{prescribed})')
    ylabel('Force PSD (dB/Hz)')
    xlim([0 10]) 
    ylim([-80 20])
%     title(['CPS ring down, simulated: m = ',num2str(M),' kg, f_{natural} = ',num2str(f_nat),' Hz, \zeta = ',num2str(zeta)]);
    hold off
