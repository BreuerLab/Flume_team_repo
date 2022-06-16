function Out = Gromit_conv(dat,Gbias)

% Calibration SI-165-15 CalDate=5/24/2019

Gromit_Cal = [-0.358090000000000,-0.0161800000000000,0.338480000000000,-17.4630700000000,0.890470000000000,16.8772400000000;0.155930000000000,21.5597800000000,-0.785970000000000,-10.0453800000000,0.430010000000000,-9.68643000000000;29.8995100000000,0.375370000000000,29.6925600000000,0.457840000000000,30.5462600000000,0.601450000000000;0.00323000000000000,0.492950000000000,-1.03800000000000,-0.243410000000000,1.06263000000000,-0.203400000000000;1.21724000000000,0.0112800000000000,-0.589940000000000,0.391800000000000,-0.631210000000000,-0.395370000000000;-0.00657000000000000,-0.716070000000000,0.0176200000000000,-0.663460000000000,-0.0425100000000000,-0.644210000000000];

dat = dat-repmat(Gbias,numel(dat(:,1)),1);
Out = dat*Gromit_Cal';