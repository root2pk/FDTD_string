%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PBMMI SIX STRING FDTD ASSIGNMENT 
%%%%% 
%%%%%
%%%%%
%%%%% FUNCTION THAT CALLS THE FDTD SINGLE STRING FUNCTION FOR SIX DIFFERENT STRINGS
%%%%% OF VARIOUS PARAMETERS, AND PRODUCES A CHORD
%%%%% THE STRINGS ARE OF A STEEL STRING GUITAR TUNED EADGBE.
%%%%% 
%%%%% PARAMETERS : 
%%%%% cho : 2X6 array with the first row containing frequency data in Hz and the
%%%%%       second row containing string lengths in metres
%%%%% rad : 1X6 array with the string gauges in inches
%%%%% type : 'up' or 'down' sets the strumming direction
%%%%% gap : Time gap between consecutive plucks in seconds
%%%%% 
%%%%% References :
%%%%% https://pages.mtu.edu/~suits/notefreqs.html
%%%%% https://www.dawsons.co.uk/blog/acoustic-guitar-strings-guide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = Chord(cho,rad,type,gap)

%% Options
opts.plot_on = false;
opts.useforloop = true;
opts.add_stiffness = true;
opts.input_type = 'struck';
opts.output_type = 'displacement';
opts.bctype = 'clamped';
%opts.bctype = 'simply_supported';


rho = 7850;                            % String density (kg/m^3)
exc_interval = gap;                    % Excitation interval (s)
Tf = 4;                                % Time of simulation (s)                 
SR = 44100;                            % Sample rate (samples/s)

%% Arrays for multiple strings

% List of string frequencies for chords (Hz)
freq_list = cho(1,:);

% List of string radii (inches)
r_list = rad;

% Converting inches to metres
r_list = r_list.*0.00254;

% Length list because string length varies while playing chords and holding
% frets (m)
L_list = cho(2,:);

% Tension list obtained using the formula ( f = (1/2L)*sqrt(T/M) ) (N)
T_list = 4*pi*rho*(L_list.^2).*(freq_list.^2).*(r_list.^2);

% Excitation time list
exc_list = [0:exc_interval:5*exc_interval];

% Excitation coordinate list
xi_list = 0.1 + 0.1*rand(1,6);

% Peak Amplitue list
famp_list = 5 + rand(1,6);

% Duration of Excitation list
dur_list = 5*rand(1,6)*0.001;

% Coordinate of Output List
xo_list = 0.6 + 0.2*rand(1,6);


%% Output Signal

% Output Signal Length

output = zeros(1,(5*exc_interval+Tf)*SR+1);  % Length of output signal in samples

parfor i = 1:6
    
    % Re initialising structs for parfor
    phys_param = struct();
    sim_param = struct();
    
    phys_param.E = 2e11;                   % Young's modulus (Pa)    
    phys_param.T60 = 5;                    % T60 (s)
    phys_param.rho = 7850;                 % density (kg/m^3)
   
    phys_param.L = L_list(i);              % length (m)
    phys_param.T = T_list(i);              % tension (N)
    phys_param.r = r_list(i);              % string radius (m)
    
    sim_param.SR = 44100;                  % sample rate (Hz)
    sim_param.Tf = 4;                      % duration of simulation (s)
    sim_param.exc_st = exc_list(i);        % start time of excitation (s)
    sim_param.xi = xi_list(i);             % coordinate of excitation (normalised, 0-1)
    sim_param.famp = famp_list(i);         % peak amplitude of excitation (N)
    sim_param.dur = dur_list(i);           % duration of excitation (s)
    sim_param.xo = xo_list(i);             % coordinate of output (normalised, 0-1)
    y(i,:) = String_FDTD(opts,phys_param,sim_param);  
    
      
end

exc_interval = exc_interval*SR;                % Excitation interval in samples

% Adding up signals into one signal
for i = 1:6
    
    start_pos = exc_interval*(i-1)+1;
    end_pos = start_pos + Tf*SR - 1;
    % Downstroke
    if strcmp(type,'down')
        output(start_pos:end_pos) = output(start_pos:end_pos) + y(i,:);
    % Upstroke
    elseif strcmp(type,'up')       
        output(start_pos:end_pos) = output(start_pos:end_pos) + y(7-i,:);
    % Block Chord
    else
        output(1:Tf*SR) = output(1:Tf*SR) + y(i,:);
    end
    
end
   


%% Play sound
soundsc(output,SR);