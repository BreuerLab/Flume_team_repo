% Calculate dataset efficiency
% Eric Handy, Jun 2022

clear;

main_dir = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Flume_team_repo\DAQandMotorControl');
cd(main_dir);

% Folderpath where files are stored

folderpath = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20220619_TandemSunday_AlphaSweep_APHPhase_A2E_a33_a67\');
% folderpath = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20220617_TandemFriday_AlphaSweep_PHPhase_A2E_a15\');

% Sweeping parameters

p3 = [60    65    70    75    80]; lp3 = length(p3);
ph = [-180  -120   -60     0    60   120   180]; lph = length(ph);
h3 = [0.5500    0.7000    0.8500    1.0000    1.1500    1.3000]; lh3 = length(h3);

% Initialize variables

Eff_sys = NaN(lp3,lph,lh3);
Eff_sys_corr = Eff_sys;

CPP2 = Eff_sys;
CPH2 = Eff_sys;
CPP3 = Eff_sys;
CPH3 = Eff_sys;

Eff_phys_2 = Eff_sys;
Eff_phys_3 = Eff_sys;

Eff_corr_2 = Eff_sys;
Eff_corr_3 = Eff_sys;

Yp = Eff_sys;

U_flow = Eff_sys;
U_2prime = Eff_sys;
U_3prime = Eff_sys;

beta2 = Eff_sys;
beta3 = Eff_sys;
alphaT4_2 = Eff_sys;
alphaT4_3 = Eff_sys;

CD2_norm = Eff_sys;
CD2_norm_p = Eff_sys;
CD3_norm = Eff_sys;
CD3_norm_p = Eff_sys;

% Calculation loop

tic;
for ii = 1:lp3
    for jj = 1:lph
        for kk = 1:lh3
            
            filename = ['20220619_TandemFoil_APHPhaseSweep_A2E_alpha=0.679_p3=',...
                num2str(p3(ii)), '_h3=', num2str(h3(kk),3), 'c_phase=', num2str(ph(jj)), '.mat'];
            
%             filename = ['20220617_TandemFoil_PHPhaseSweep_A2E_p3='...
%                 num2str(p3(ii)), '_h3=', num2str(h3(kk),3), 'c_phase=', num2str(ph(jj)), '.mat'];
            
            filepath = fullfile(folderpath,filename);
            load(filepath);
            
            out(:,5) = deg2rad(Prof_out_angle(:,5)); % for data taken on 20220619
            
            [kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, 1000, 3, 6*0.0762); % extract parameters and
            res = calculate_forces(par, kin, out);                                                           % perform calculations
            
            % store calculated values
            
            CPP2(ii,jj,kk) = mean(res.CPP2);
            CPH2(ii,jj,kk) = mean(res.CPH2);
            CPP3(ii,jj,kk) = mean(res.CPP3);
            CPH3(ii,jj,kk) = mean(res.CPH3);
            
            Eff_sys(ii,jj,kk) = res.Eff_sys;
            Eff_sys_corr(ii,jj,kk) = res.Eff_sys_corr;
            
            Eff_phys_2(ii,jj,kk) = res.Eff_2;
            Eff_phys_3(ii,jj,kk) = res.Eff_3;
            
            Eff_corr_2(ii,jj,kk) = res.Eff_2prime;
            Eff_corr_3(ii,jj,kk) = res.Eff_3prime;
            
            Yp(ii,jj,kk) = res.Yp;
            
            U_flow(ii,jj,kk) = par.U;
            U_2prime(ii,jj,kk) = res.U_2prime;
            U_3prime(ii,jj,kk) = res.U_3prime;
            
            beta2(ii,jj,kk) = res.beta2;
            beta3(ii,jj,kk) = res.beta3;
            alphaT4_2(ii,jj,kk) = par.alphaT4_2;
            alphaT4_3(ii,jj,kk) = par.alphaT4_3;
            
            CD2_norm(ii,jj,kk) = res.CD2_norm;
            CD2_norm_p(ii,jj,kk) = res.CD2_norm_p;
            CD3_norm(ii,jj,kk) = res.CD3_norm;
            CD3_norm_p(ii,jj,kk) = res.CD3_norm_p;
            
        end
    end
end
toc;

save('20220726_TandemFoil_calcWake_efficiency_A2E_a68_PHPh.mat',...
    'foiltype', 'p3', 'h3', 'ph', 'alphaT4', 'CPP2', 'CPH2', 'CPP3', 'CPH3', 'Eff_sys', 'Eff_sys_corr',...
    'Eff_phys_2', 'Eff_phys_3', 'Eff_corr_2', 'Eff_corr_3', 'Yp', 'U_flow', 'U_2prime', 'U_3prime', 'beta2', 'beta3',...
    'alphaT4_2', 'alphaT4_3', 'CD2_norm', 'CD2_norm_p', 'CD3_norm', 'CD3_norm_p');

