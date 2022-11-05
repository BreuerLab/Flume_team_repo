%% Testing code, selecting points in a figure

clear;

piv_freq = 14.9316; % frequency of piv data
frame = 23; % piv frames per cycle

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=16_p3=75_h3=0.7_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU_01\Export\TandemFoil_aT4=0.16_p3=75_h3=0.7c_ph=-120.mat');

%% coordinate selection

for ii = 1:frame
    
    figure('Name', 'Foil Coordinates', 'WindowState', 'maximized');
    contour(compiled_data(ii).x, compiled_data(ii).y, compiled_data(ii).isValid);
    colormap('bone');
    fprintf(['\n\nSelect approximate foil center coordinates for the LEADING and TRAILING foil in that order.\n\n When done with the current frame press [Enter]\n\nCurrent frame: ',num2str(ii),'\n']);
    leading = drawpoint('Label','Leading Foil');
    trailing = drawpoint('Label','Trailing Foil');
    coords_le(ii,:) = leading.Position;
    coords_tr(ii,:) = trailing.Position;
    close 'Foil Coordinates';
    
end