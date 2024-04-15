% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % %     CHORD GENERATOR WITH CUSTOMISATION
% % %     RUTHU PREM KUMAR 
% % %     MARCH 2020
% % %     
% % %     
% % %     THIS PROGRAM CALLS A FUNCTION WHICH USES THE SINGLE STRING FDTD
% % %     FUNCTION TO GENERATE CHORDS
% % %     THE PARAMETERS THAT CAN BE CHANGED ARE
% % %     
% % %     THE CHORD ITSELF
% % %     THE GAUGE OF THE STEEL STRINGS
% % %     THE METHOD OF STRIKING (DOWN OR UP)
% % %     THE GAP BETWEEN CONSECUTIVE STRIKES
% % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear all; close all;

%% String length when fret is held (in metres)

f0 = 0.6477;       % open
f1 = 0.6113;       % 1st fret
f2 = 0.5771;       % 2nd fret      
f3 = 0.5447;       % 3rd fret
f4 = 0.5142;       % 4th fret
f5 = 0.4854;       % 5th fret

%% Frequencies (In Hz)
% This forms a 2X6 matrix with the first row providing frequency data and
% the second row providing string length data

Am_chord = [82.41,110,164.81,220,261.63,329.63;f0,f0,f2,f2,f1,f0];
C_chord = [82.41,130.81,164.81,196,261.63,329.63;f0,f3,f2,f0,f1,f0];
G_chord = [98,123.47,146.83,196,293.66,392;f3,f2,f0,f0,f3,f3];
Em_chord = [82.41,123.47,164.81,196,246.94,329.63;f0,f2,f2,f0,f0,f0];
F_chord = [87.31,130.81,174.61,220,261.63,349.23;f1,f3,f3,f2,f1,f1];
Dm_chord = [87.31,110,146.83,220,293.66,349.23;f1,f0,f0,f2,f3,f1];
D_chord = [92.50,110,146.83,220,293.66,369.99;f2,f0,f0,f2,f3,f2];
Bm_chord = [92.50,123.47,185,246.94,293.66,369.99;f2,f2,f4,f4,f3,f2];
Bflat_chord = [87.31,116.54,174.61,233.08,293.66,349.23;f1,f1,f3,f3,f3,f1];
Gm_chord = [98,146.83,196,233.08,293.66,392;f3,f5,f5,f3,f3,f3];

%% String Gauges (in inches)
% From lightest to heaviest

Extra_Light = [0.010,0.014,0.023,0.030,0.039,0.047];
Custom_Light = [0.011,0.015,0.023,0.032,0.042,0.052];
Light = [0.012,0.016,0.025,0.032,0.042,0.054]; 
Medium = [0.013,0.017,0.026,0.035,0.045,0.056]; 
Heavy = [0.014,0.018,0.027,0.039,0.049,0.059]; 

% From heaviest to lightest
Extra_Light = flip(Extra_Light);
Custom_Light = flip(Custom_Light);
Light = flip(Light);
Medium = flip(Medium);
Heavy = flip(Heavy);
%%%%%%%%%%%%%%%% -------------------%%%%%%%%%%%%%%%%%%%%

%% Call chord function
% Format - Chord(<chord name>,<gauge name>,'up'/'down',strike gap(s))
% If neither up nor down is specified, the function automatically plays all
% the strings together as a block chord
Chord(Am_chord,Heavy,'down',0.1);
