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
Eff_sys_corr = NaN(lp3,lph,lh3);
Eff_phys_2 = NaN(lp3,lph,lh3);
Eff_phys_3 = NaN(lp3,lph,lh3);
Eff_corr_2 = NaN(lp3,lph,lh3);
Eff_corr_3 = NaN(lp3,lph,lh3);

alphaT4_2 = NaN(lp3,lph,lh3);
alphaT4_3 = NaN(lp3,lph,lh3);
CD2_norm = NaN(lp3,lph,lh3);
CD3_norm = NaN(lp3,lph,lh3);
CD2_norm_p = NaN(lp3,lph,lh3);
CD3_norm_p = NaN(lp3,lph,lh3);
Yp = NaN(lp3,lph,lh3);
U_flow = NaN(lp3,lph,lh3);
U_2prime = NaN(lp3,lph,lh3);
U_3prime = NaN(lp3,lph,lh3);

% Calculation loop

tic;
for ii = 1:lp3
    for jj = 1:lph
        for kk = 1:lh3
            
            filename = ['20220619_TandemFoil_APHPhaseSweep_A2E_alpha=0.33_p3=',...
                num2str(p3(ii)), '_h3=', num2str(h3(kk),3), 'c_phase=', num2str(ph(jj)), '.mat'];
            
%             filename = ['20220617_TandemFoil_PHPhaseSweep_A2E_p3='...
%                 num2str(p3(ii)), '_h3=', num2str(h3(kk),3), 'c_phase=', num2str(ph(jj)), '.mat'];
            
            filepath = fullfile(folderpath,filename);
            load(filepath);
            
            out(:,5) = deg2rad(Prof_out_angle(:,5)); % for data taken on 20220619
            
            [kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, 1000, 3, 6*0.0762); % extract parameters and
            res = calculate_forces(par, kin, out);                                                               % perform calculations
            
            % store calculated values
            
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

save('20220706_BarnsleyWellicome_TandemFoil_efficiency_A2E_a33_PHPh.mat',...
    'foiltype', 'p3', 'h3', 'ph', 'alphaT4', 'Eff_sys', 'Eff_sys_corr', 'Eff_phys_2', 'Eff_phys_3', 'Eff_corr_2', 'Eff_corr_3', 'U_flow', 'alphaT4_2', 'alphaT4_3', 'CD2_norm', 'CD2_norm_p', 'CD3_norm', 'CD3_norm_p', 'Yp', 'U_2prime', 'U_3prime');

