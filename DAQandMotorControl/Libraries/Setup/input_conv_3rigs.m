function Prof_out = input_conv_3rigs(profs, freq, heave1, heave2, heave3,pitch_bias)
% Input is [Profiles (Shawn Wallace Gromit),frequency, heave amplitdues ]
% Input needs to be Nx6 (technically)
% Each rig must have a pitch column [deg] and a heave column [meters]

% Prof_out = input_conv_3rigs(Prof_out_angle, freq1, heave1, heave2, heave3);

if size(profs,2)~=6
    error('Input must be Nx6!')
end

%% Heave 1 and 2 (Shawn and Gromit)

% out(:,2) = -profs(:,2)./0.03*1.0590; % why is this here?
% out(:,4) = -profs(:,4)./0.03/0.9891;

%   f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2
%   Coefficients (with 95% confidence bounds):

       p00 =  0.9543*1.975/2;   %  (  0.9437  ,  0.965   )
       p10 = -0.002933;         %  ( -0.02833 ,  0.02246 )
       p01 = -0.03169;          %  ( -0.05104 , -0.01234 )
       p20 = -0.02473;          %  ( -0.04419 , -0.005266)
       p11 = -9.101e-05;        %  ( -0.01196 ,  0.01177 )
       p02 =  0.01964;          %  (  0.008342,  0.03094 )
       
       x = freq;
       y = heave1;
       
    Prof_out(:,2) = -profs(:,2)./0.03./(p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2); % Heave 1 (Shawn)
    
%   f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2
%   Coefficients (with 95% confidence bounds):

       p00 =  0.996;    %  (  0.9921  ,  0.9999  )
       p10 =  0.002508; %  ( -0.006756,  0.01177 )
       p01 = -0.135;    %  ( -0.2276  , -0.04232 )
       p20 = -0.01149;  %  ( -0.01859 , -0.004386)
       p11 = -0.03915;  %  ( -0.09595 ,  0.01765 )
       p02 =  1.232;    %  (  0.5224  ,  1.942   )
       
       x = freq;
       y = heave2;
       
    Prof_out(:,4) = -profs(:,4)./0.03/(p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2); % Heave 2 (Gromit)

%% Heave 3 (Wallace)

%   f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y
%   Coefficients (with 95% confidence bounds):

%        p00 =  0.9793;  % (  0.9755  ,  0.9831  )
%        p10 =  0.04503; % (  0.03297 ,  0.0571  )
%        p01 = -0.01027; % ( -0.05548 ,  0.03494 )
%        p20 = -0.0365;  % ( -0.04597 , -0.02704 )
%        p11 = -0.01869; % ( -0.08684 ,  0.04946 )

%    f(x,y) = p00 + p10*x + p01*y
%    Coefficients (with 95% confidence bounds):

       p00 =  0.991;    %  (  0.0001948,  0.0006857 )
       p10 = -0.01095;  %  ( -0.001101 , -0.000552  )
       p01 =  0.04966;  %  (  1.009    ,  1.013     )
       
       x = freq;
       y = heave3;
       
    Prof_out(:,6) = -profs(:,6)./0.03/(p00+p10*x + p01*y); % Heave 3 (Wallace)

%% Pitch 1 and 2 (Shawn and Gromit)

Prof_out(:,1) = -profs(:,1).*5./(360)*2*45/56*5 + pitch_bias(1); % Pitch 1 (Shawn)
% Prof_out(:,3) = -profs(:,3).*5./(360)*2*60/74 + pitch_offset2;%/1.0126
% below is for Parker Motion motor
% Prof_out(:,5) = -profs(:,5)*5.*10/30000*12800/360*1.00+pitch_bias(3); % 5:1 gear *( 10 volts /30000 steps) * (12800 steps/ revolution) * calibration gain 4/4/2016
Prof_out(:,5) = -profs(:,5)*5.*0.022222+pitch_bias(3); % 5:1 gear *(1 Volt/ 184.613 deg) for new Teknic servo motor
% ^^^^^ Supposedly this is sent out to Wallace?

%% Pitch 3 (Servo motor)

% Prof_out(:,3) = -profs(:,3)*1.*10/2000*1000/360*1.00*90/57.42+pitch_offset2;
Prof_out(:,3) = -profs(:,3).*5./(360)*2*65/81*5 + pitch_bias(2);
% ^^^^^ and this is sent out to Gromit?

if abs(Prof_out)>10*ones(numel(Prof_out(:,1)),6)
    error('Voltage output too high! Check conversion')
end


end