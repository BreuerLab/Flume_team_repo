%% run ldv test - 20220930

[foil, ~, ~] = foils_database('C1');
chord = foil.chord;

loops = 350;

total_predicted_time = loops*((1/0.6492)*(40 + 2*3) + 20)/60/60

fprintf('press any key to run... \n')

pause

fprintf('running...\n')

for i = 1:loops

    [flume, out, dat, Prof_out_angle, Prof_out, last_out, freq, pitch2, heave2, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
    = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,...
    'C1', 0.6492, 40, 1, 40, 1, -180, -90, 40, 3, 0, 0, 0, 4);

    file = ['20220930_cylinder_ldv_test_', num2str(i), '.mat'];
    save(file);
    
    fprintf(sprintf('Loop: %i \n', i));

    pause(20);

end