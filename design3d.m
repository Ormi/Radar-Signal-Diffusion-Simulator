%%
% @author xormos00
% @date March 2017
% @title Radar Signal Diffusion Simulator
% @Bachelor Thesis
% VUT FIT
%
% @dependecies
% return_signal_freq.m
% radar_equation.m
% parse_csv_file.m
% make_plot.m
%
% @using HGtransform
% @literature Matlab
%
% Version 1.0

% JSONLab
% oprav hlupe nazvy premmennych a funkcii, rovnako ako aj funkcii
% posun abs
% posun phi init
% chyba v generovani
%
%


% @TODO 
% RCS pre kazdy bod jedinecne
% Spectogram
% Definovat poziciu bodov ako trojicu a nasobit medzi

% @QUESTIONS
% -Aky RCS urcit pre jednotlive body?

% @IDEAS
% GUI like RadarWaveformAnalyzer?
%   https://www.mathworks.com/help/phased/ug/radar-waveform-analyzer-app.html

clear all;

%%
% Setting static variables for simulation

F_TRANS = 24.125e+8;           % Transmiting frequency GHz
SPEED_OF_OBJECT = 30;       % Speed of object
ANTENNA_GAIN = 13;          % dB 

OBJECT_X_POS = 22;          %   __
OBJECT_Y_POS = 56;          % _/  \__
OBJECT_Z_POS = 0;           %

RADAR_X_POS = 15;
RADAR_Y_POS = 20;
RADAR_Z_POS = 15;
RADAR_AIM_X_POS = 15;
RADAR_AIM_Y_POS = 40;
RADAR_AIM_Z_POS = 0;


NUM_OF_STEPS = 1e4;        % Simulation run repetitions (frequency)
                            % for best results us 10kHz

SHOW_CUBE_SIMULATION = 0;   % Put 1 fow show whole cube Simulation
                            % Time demanding action!! Graphical interface

%%
% Setting up size of enviroment(cube) for simulation
ah = axes;
set(ah,'XLim',[0 30],'YLim',[0 70],...
    'ZLim',[0 30]);
hold on;
view(3);

%%
% Reading and procesing files for antenna system diagram loss
antenna_hori_data = parse_csv_file('KMC4_antena_char_hori.csv');
antenna_vert_data = parse_csv_file('KMC4_antena_char_vert.csv');

%%
% Enviroment
% Moving points for creating object shape.. 
%  o o
% o   o o

% Generating moving object - car
directions = [0 0 2; 0 -2 0; 0 2 2; 0 4 0; 0 -4 0; 4 0 2; 4 -2 0; 4 2 2; 4 4 0; 4 -4 0; 2 -4 0; 2 4 0];
for num=1:12
    hpoint(num) = line('XData', OBJECT_X_POS+directions(num,1),'YData', OBJECT_Y_POS+directions(num,2),'ZData', OBJECT_Z_POS+directions(num,3),'Color','black','Marker',...
       'o','MarkerSize',2,'MarkerFaceColor','black');
end

% Static point - radar
spoint = line('XData', RADAR_X_POS,'YData', RADAR_Y_POS,'ZData', RADAR_Z_POS,'Color','red','Marker',...
    'x','MarkerSize',10,'MarkerFaceColor','red');
spoint_aim = line('XData', RADAR_AIM_X_POS,'YData', RADAR_AIM_Y_POS,'ZData', RADAR_AIM_Z_POS,'Color','red','Marker',...
    'x','MarkerSize',10,'MarkerFaceColor','red');

% Library class for object movement
ht = hgtransform('parent',ah);
set(hpoint,'Parent',ht);

%%
% SIMULATION
Fs = NUM_OF_STEPS;
dt = 1/Fs;
t = 0:dt:1;
for i=1:length(t);
    % Object movement in enviroment settinhs
    tx = 0;
    ty = -(25/NUM_OF_STEPS)*i*2;
    tz = 0;
    
    % Transform enviroment(move object) in every step
    trans = makehgtform('translate',[tx ty tz]);
    set(ht,'Matrix',trans);

    % Show transofrmation in every step on/off
    if (SHOW_CUBE_SIMULATION)
        pause(0.1);
    end
    
    % Get values from static points
    s1x = get(spoint,'XData'); s1y = get(spoint,'YData'); s1z = get(spoint,'ZData');
    s2x = get(spoint_aim,'XData'); s2y = get(spoint_aim,'YData'); s2z = get(spoint_aim,'ZData');    
    %h1x = get(hpoint(1),'XData'); h1y = get(hpoint(1),'YData'); h1z = get(hpoint(1),'ZData');

    % Get valus from dynamic points
    % For better quicker calculations in separately cycle
    for pnt=1:12
        h1(pnt) = get(hpoint(pnt),'XData'); h1y = get(hpoint(pnt),'YData'); h1z = get(hpoint(pnt),'ZData');     
    end

    %%
    % In every step of simulation we have to calculate data for every moving point
    for pnt=1:12          
        h1x = h1(pnt);   
        
        % Calculation of distance
        distance(i) = abs(norm([s1x, s1y, s1z] - [h1x+tx, h1y+ty, h1z+tz]));
        
        % Point_start to Point_actual_position % start_point - actual_point
        heading_p2p = [h1x - (h1x + tx), h1y - (h1y + ty), h1z - (h1z + tz)];
        
        % Radar to Point_actual_position % start_point - actual_point
        heading_r2p = [s1x - (h1x + tx), s1y - (h1y + ty), s1z - (h1z + tz)];
        heading_p2r = [(h1x + tx) - s1x, (h1y + ty) - s1y, (h1z + tz) - s1z];    

        % Radar to Radar looking point % start_point - actual_point
        heading_r2r = [s1x - s2x, s1y - s2y, s1z - s2z];

        %angle1 = [abs(heading_r2p(1)), abs(heading_r2p(2)), abs(heading_r2p(3))]; 
        vector_p2r = [heading_p2r(1), heading_p2r(2), heading_p2r(3)];     
        vector_p2p = [heading_p2p(1), heading_p2p(2), heading_p2p(3)]; 
        % Area angle
        angle_area = atan2(norm(cross(vector_p2r,vector_p2p)),dot(vector_p2r,vector_p2p));
      
        % Vertical calculation of angle between radar vector and radar2object vector
        % x axis is neglected for 3D to 2D transformation 
        vector_r2p_temp = [0, (heading_r2p(2)), (heading_r2p(3))];     
        vector_r2r_temp = [0, (heading_r2r(2)), (heading_r2r(3))]; 
        angle_vert(i) = rad2deg(atan2(norm(cross(vector_r2p_temp,vector_r2r_temp)),dot(vector_r2p_temp,vector_r2r_temp)));  
        
        % Horizontal calculation of angle between radar vector and radar2object vector
        % z axis is neglected for 3D to 2D transformation
        vector_r2p_temp = [(heading_r2p(1)), (heading_r2p(2)), 0];
        vector_r2r_temp = [(heading_r2r(1)), (heading_r2r(2)), 0]; 
        angle_hori(i) = rad2deg(atan2(norm(cross(vector_r2r_temp,vector_r2p_temp)),dot(vector_r2r_temp,vector_r2p_temp)));    

        % Angle correction
        % In .csv input file we have angles from 0-180 degress, we have to make a correction
        ang1 = 9e+4 + angle_hori(i) * 1e+4;
        ang2 = 9e+4 + angle_vert(i) * 1e+4;

        % Get loss from antenna system diagram
        antenna_loss_hori(i) = antenna_hori_data(logical(ang1));
        antenna_loss_vert(i) = antenna_vert_data(logical(ang2));

        antenna_loss = antenna_loss_vert(i) * antenna_loss_hori(i);

        % #Extern function for calculating receiving frequency
        % @input is speed of object and transmiting frequency
        F_receiv(pnt,i) = return_signal_freq(SPEED_OF_OBJECT,F_TRANS, angle_area);

        % #Extern function for calculating receiving power
        % @input is speed of object and transmiting frequency
        P_receiv(pnt,i) = radar_equation(F_TRANS, 1, F_TRANS, antenna_loss, distance(i));
    end

    if mod(i, NUM_OF_STEPS/100) == 0
        disp('Simulation completition status: [in %]');
        disp(100/(NUM_OF_STEPS/i));
    end


end 

%%
% Generating signal
Fs = NUM_OF_STEPS;
dt = 1/Fs;
t = 0:dt:1;

for pnt=1:12
    phi_t1 = 0;    
	for n=1:length(t)
	        phi_t2 = 2*pi*F_receiv(pnt,n)*dt;
	        var = sqrt(P_receiv(pnt,n));
	        x(pnt,n) = var*exp(1j*(phi_t1+phi_t2)); % complex signal
	        phi_t1 = phi_t2;
	end

    disp('Calculations completition status: [in %]');
    disp((100/12)*pnt);    
end

%%
% Link frequency characteristic of all points together
for n=1:NUM_OF_STEPS
    x(1, n) = x(1, n) + x(2, n) + x(3, n) + x(4, n) + x(5, n) + x(6, n) + x(7, n) + x(8, n) + x(9, n) + x(10, n) + x(11, n) + x(12, n);
end

%%
make_plot(x(1,:),NUM_OF_STEPS);
