%% Belt Traverse Testing

foiltype = 'A3E';
[foil, ~, ~] = foils_database(foiltype);

freq = 0.6;
pitch2 = 75; % traverse (downstream)
pitch3 = 75; % wallace (upstream)
heave2 = 1*foil.chord;
heave3 = 1*foil.chord;
phase13 = 0;
phi = -90;
num_cyc = 40;
transientcycs = 5;
constantamplitude = 0;
offset = 0;

[flume, out, dat, Prof_out_angle, Prof_out,last_out, freq, pitch2, heave2, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
    = run_Motors(dq,last_out,bias,foiltype, freq, pitch2, heave2, pitch3, heave3, phase13, phi,...
    num_cyc, transientcycs, constantamplitude, offset);

save('test_20230216_new_traverse');