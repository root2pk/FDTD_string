%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PBMMI SIX STRING FDTD ASSIGNMENT 
%%%%% 
%%%%%
%%%%%
%%%%% PROGRAM THAT CALLS THE FDTD SINGLE STRING FUNCTION FOR SIX DIFFERENT STRINGS
%%%%% OF VARIOUS PARAMETERS.
%%%%% THE STRINGS ARE OF A STEEL STRING GUITAR TUNED EADGBE.
%%%%%
%%%%% THE OUTPUT IS A STEREO OUTPUT WITH EACH STRING AT DIFFERENT PAN
%%%%% LEVELS
%%%%%
%%%%% References : 
%%%%% https://pages.mtu.edu/~suits/notefreqs.html
%%%%% https://www.dawsons.co.uk/blog/acoustic-guitar-strings-guide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%% Options
opts.plot_on = false;
opts.useforloop = true;
opts.add_stiffness = false;
opts.input_type = 'plucked';
opts.output_type = 'displacement';
%opts.bctype = 'clamped';
opts.bctype = 'simply_supported';

%% Fixed Physical Parameters 

phys_param.E = 2e11;                   % Young's modulus (Pa)
phys_param.L = 0.6477;                 % length (m)
phys_param.T60 = 5;                    % T60 (s)
phys_param.rho = 7850;                 % density (kg/m^3)
exc_interval = 0.5;                    % Excitation interval (s)

%% Arrays for multiple strings

% List of string frequencies (E2,A2,D3,G3,B3,E4)
freq_list = [82.41,110,146.83,196,246.94,329.63];

% List of string radii (Edit based on gauge values in inches)
r_list = [0.054,0.042,0.032,0.025,0.016,0.012];

% Converting inches to metres
r_list = r_list.*0.00254;

% Tension list obtained using the formula ( f = (1/2L)*sqrt(T/M) )
T_list = 4*pi*phys_param.rho*(phys_param.L^2)*(freq_list.^2).*(r_list.^2);

% Excitation time list
exc_list = [0:exc_interval:5*exc_interval];

% Excitation coordinate list
xi_list = [0.7,0.8,0.85,0.75,0.8,0.77];

% Peak Amplitue list
famp_list = [5,5,5,5,5,5];

% Duration of Excitation list
dur_list = [0.001,0.002,0.001,0.003,0.002,0.001];

% Coordinate of Output List
xo_list = [0.1,0.1,0.1,0.1,0.1,0.1];

% Stereo list
stereo_list = linspace(-1,1,6);

%% Fixed Simulation parameters 
sim_param.SR = 44100;                  % sample rate (Hz)
sim_param.Tf = 4;                      % duration of simulation (s)


%% Output Signal
% Output Signal Length

output_L = zeros((5*exc_interval + sim_param.Tf)*sim_param.SR+1,1);  % Length of left output signal in samples
output_R = output_L;                                                 % Length of right output signal in samples
exc_interval = exc_interval*sim_param.SR;                            % Excitation interval in samples
for i = 1:6
    %% Variable physical parameters 
    phys_param.T = T_list(i);          % tension (N)
    phys_param.r = r_list(i);          % string radius (m)
    sim_param.exc_st = exc_list(i);    % start time of excitation (s)
    sim_param.xi = xi_list(i);         % coordinate of excitation (normalised, 0-1)
    sim_param.famp = famp_list(i);     % peak amplitude of excitation (N)
    sim_param.dur = dur_list(i);       % duration of excitation (s)
    sim_param.xo = xo_list(i);         % coordinate of output (normalised, 0-1)
    y = String_FDTD(opts,phys_param,sim_param); 
    
    %% Stereo properties
    
    % Contribution to each channel based on pan level
    coeff1 = abs(1 - stereo_list(i));  
    coeff2 = abs(-1 - stereo_list(i));    
    y_L = coeff1*y;
    y_R = coeff2*y;  
    
    % Adding all the signals into two channels
    start_pos = exc_interval*(i-1)+1;
    end_pos = start_pos + sim_param.Tf*sim_param.SR - 1;
    output_L(start_pos:end_pos) = output_L(start_pos:end_pos) + y_L;
    output_R(start_pos:end_pos) = output_R(start_pos:end_pos) + y_R;
    
end

% Converting two channels to one stereo signal
output = [output_L,output_R];

%% Play sound
soundsc(output,sim_param.SR);
