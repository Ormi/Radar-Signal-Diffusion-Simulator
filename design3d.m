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
% KMC4_antena_char_hori.csv
% KMC4_antena_char_vert.csv
% /jsonlab
% scene.json
% model_XXX.json

% @using HGtransform
% @literature Matlab
%
% Version 1.0

% @TODO 
% 

% @DONE
% rozdelenie suborov, vytvorenie viacerych modelov.

% @QUESTIONS
% Nemam vela kapitol ?
% Dokoncit 2 kanaly.
% .json suboru musia byt 2, nejdu len tak lahko do seba includovat.
% finalna dokumentacia, chcete ju este citat, alebo ju mozme volne odovzdat?

clear all;
addpath(genpath('./jsonlab'))
savepath

%%
% Setting static variables for simulation
data_scene=loadjson('scene.json');
data_model=loadjson('model_truck.json');

% Transmiting frequency GHz
F_TRANS = str2num(data_scene.global.Transceiving_frequency);
% Speed of object
SPEED_OF_OBJECT = str2num(data_scene.model.Speed_of_object);
% dB 
ANTENNA_GAIN = str2num(data_scene.global.Antenna_gain);

OBJECT_POS = [str2num(data_scene.model.position.x) str2num(data_scene.model.position.y) str2num(data_scene.model.position.z)];

RADAR_POS = [str2num(data_scene.radar.position.x) str2num(data_scene.radar.position.y) str2num(data_scene.radar.position.z)];

RADAR_AIM = [str2num(data_scene.radar.heading.x) str2num(data_scene.radar.heading.y) str2num(data_scene.radar.heading.z)];

ENVIROMENT_SIZE = [str2num(data_scene.enviroment.x) str2num(data_scene.enviroment.y) str2num(data_scene.enviroment.z)];

% Simulation run repetitions (frequency)
% for best results us 10kHz
NUM_OF_STEPS = str2num(data_scene.global.Number_of_steps);

RCS = str2num(data_scene.model.RCS);

% Put 1 fow show whole cube Simulation
% Time demanding action!! Graphical interface
Simulation_view = (data_scene.global.Show_simulation_movement);
if strcmp(Simulation_view,'on')
    SHOW_CUBE_SIMULATION = 1;
else
    SHOW_CUBE_SIMULATION = 0;
end

Simulation_setup = (data_scene.global.Test_simulation_setup);
if strcmp(Simulation_setup,'on')
    SHOW_SETUP = 1;
else
    SHOW_SETUP = 0;
end

Radar_rotation = (data_scene.global.Right_angle_rotation_z_axis);
if strcmp(Radar_rotation,'on')
    SWITCH_AXIS = 1;
else
    SWITCH_AXIS = 0;
end


c = physconst('Lightspeed');

%%
% Setting up size of enviroment(cube) for simulation
ah = axes;
set(ah,'XLim',[0 ENVIROMENT_SIZE(1)],'YLim',[0 ENVIROMENT_SIZE(2)],...
    'ZLim',[0 ENVIROMENT_SIZE(3)]);
hold on;
view(3);

%%
% Reading and procesing files for antenna system diagram loss
if (SWITCH_AXIS)
    antenna_vert_data = parse_csv_file('KMC4_antena_char_hori.csv');
    antenna_hori_data = parse_csv_file('KMC4_antena_char_vert.csv');
else
    antenna_hori_data = parse_csv_file('KMC4_antena_char_hori.csv');
    antenna_vert_data = parse_csv_file('KMC4_antena_char_vert.csv');
end
%%
% Enviroment
% Moving points for creating object shape.. 
%  o o
% o   o o

% Generating moving object
directions = (data_model.directions);
NUM_OF_POINTS = numel(directions)/3

for num=1:NUM_OF_POINTS
    hpoint(num) = line('XData', OBJECT_POS(1)+directions(num,1),'YData', OBJECT_POS(2)+directions(num,2),'ZData', OBJECT_POS(3)+directions(num,3),'Color','black','Marker',...
       '>','MarkerSize',2,'MarkerFaceColor','black');
end

% Static point - radar
spoint = line('XData', RADAR_POS(1),'YData', RADAR_POS(2),'ZData', RADAR_POS(3),'Color','red','Marker',...
    '*','MarkerSize',10,'MarkerFaceColor','red');
spoint_aim = line('XData', RADAR_AIM(1),'YData', RADAR_AIM(2),'ZData', RADAR_AIM(3),'Color','red','Marker',...
    'x','MarkerSize',10,'MarkerFaceColor','red');

% Library class for object movement
ht = hgtransform('parent',ah);
set(hpoint,'Parent',ht);

%%
% SIMULATION
Fs = NUM_OF_STEPS;
dt = 1/Fs;
t = 0:dt:1;

if (SHOW_SETUP)
    t = 1;
end

for i=1:length(t);
    % Object movement in enviroment settings
    tx = 0;
    ty = -((ENVIROMENT_SIZE(2)/2)/NUM_OF_STEPS)*i*1.1;
    tz = 0;
    
    % Transform enviroment(move object) in every step
    trans = makehgtform('translate',[tx ty tz]);
    set(ht,'Matrix',trans);

    % Show transofrmation in every step on/off
    if (SHOW_CUBE_SIMULATION) or (SHOW_SETUP)
        pause(0.1);
    end
    
    % Get values from static points
    s1x = get(spoint,'XData'); s1y = get(spoint,'YData'); s1z = get(spoint,'ZData');
    s2x = get(spoint_aim,'XData'); s2y = get(spoint_aim,'YData'); s2z = get(spoint_aim,'ZData');    
    %h1x = get(hpoint(1),'XData'); h1y = get(hpoint(1),'YData'); h1z = get(hpoint(1),'ZData');

    % Get valus from dynamic points
    % For better quicker calculations in separately cycle
    for point_init=1:NUM_OF_POINTS
        %h1_x(point_init) = get(hpoint(point_init),'XData'); 
        h1_x(point_init) = get(hpoint(point_init),'XData');         
        h1_y(point_init) = get(hpoint(point_init),'YData'); 
        h1_z(point_init) = get(hpoint(point_init),'ZData');     
    end

    %%
    % In every step of simulation we have to calculate data for every moving point
    for point=1:NUM_OF_POINTS         
        h1x = h1_x(point); 
        h1y = h1_y(point); 
        h1z = h1_z(point);  
         
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
        angle_area(i) = atan2(norm(cross(vector_p2r,vector_p2p)),dot(vector_p2r,vector_p2p));

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
        ang_temp1(i) = 9e+4 + (angle_hori(i) * 1e+3);
        search_angle1(i) = round(ang_temp1(i));

        % Get loss from antenna system diagram
        search_angle1(i) = 90000;
        antenna_loss_hori(i) = antenna_hori_data((search_angle1(i)));

        %disp(angle_vert);
        ang_temp2(i) = 9e+4 + (angle_vert(i) * 1e+3);   
        search_angle2(i) = round(ang_temp2(i));
        antenna_loss_vert(i) = antenna_vert_data((search_angle2(i)));

        % objekt y <= radar y
        % sme mimo pohlad radu = -70 loss
        antenna_loss(i) = (antenna_loss_vert(i) * antenna_loss_hori(i))*1e4;

        % #Extern function for calculating receiving frequency
        % @input is speed of object and transmiting frequency
        F_receiv(point,i) = return_signal_freq(SPEED_OF_OBJECT,F_TRANS, angle_area(i), c);

        % #Extern function for calculating receiving power
        % @input is speed of object and transmiting frequency
        P_receiv(point,i) = radar_equation(F_TRANS, c, RCS, antenna_loss(i), distance(i));
    end

    if mod(i, NUM_OF_STEPS/100) == 0
        disp('Simulation completition status: [in %]');
        disp(100/(NUM_OF_STEPS/i));
    end


end


if (SHOW_SETUP)
    disp('SETUP VIEW DONE');
else
    %%
    % Generating signal
    Fs = NUM_OF_STEPS;
    dt = 1/Fs;
    t = 0:dt:1;

    for point=1:NUM_OF_POINTS
        phi_t1 = 0;    
    	for n=1:length(t)
            phi_dt = 2*pi*F_receiv(point,n)*dt;%+phi
            var = sqrt(P_receiv(point,n));
            phi_t2 = phi_t1+phi_dt;
            x(point,n) = var*exp(1j*(phi_t2)); % complex signal
            phi_t1 = phi_t2;			
    	end

        disp('Calculations completition status: [in %]');
        disp((100/12)*point);    
    end

    for point=1:NUM_OF_POINTS
        phi_t1 = 0;    
        for n=1:length(t)
            phi_dt = 2*pi*F_receiv(point,n)*dt+phi_t1;
            var = sqrt(P_receiv(point,n));
            phi_t2 = phi_t1+phi_dt;
            xx(point,n) = var*exp(1j*(phi_t2)); % complex signal
            phi_t1 = phi_t2;            
        end

        disp('Calculations completition status: [in %]');
        disp((100/12)*point);    
    end


    %%
    % Link frequency characteristic of all points together
    for n=1:NUM_OF_STEPS
        for m=2:NUM_OF_POINTS
           x(1, n) = x(1, n) + x(m, n); 
           xx(1, n) = xx(1, n) + xx(m, n);
        end
        %x(1, n) = x(1, n) + x(2, n) + x(3, n);
        %x(1, n) = x(1, n) + x(2, n) + x(3, n) + x(4, n) + x(5, n) + x(6, n) + x(7, n) + x(8, n) + x(9, n) + x(10, n) + x(11, n) + x(12, n);
    end

    %%
    %make_plot(x(1,:),xx(1,:),NUM_OF_STEPS);

    %clear; clc;
    % t = 1:NUM_OF_STEPS+1;
    % x = x(1,:);
    % y = xx(1,:);
    % plot(t, real(x), 'r-',t, real(y), 'b-');

    sig = x(1, :);

    Fs=NUM_OF_STEPS;
    t = (1:length(sig))/Fs;

    SEG_LEN=128;
    OVERLAP=0; %SEG_LEN/2;
    FFT_LEN=8*SEG_LEN;

    seg_idx=1:(SEG_LEN-OVERLAP):length(sig);

    Xdb = [];

    % loop over segments and process
    for i=1:(length(seg_idx)-1)
        seg = sig(seg_idx(i):(seg_idx(i)+SEG_LEN-1));
        seg = seg - mean(seg);
        X = fftshift(fft(seg.*hamming(SEG_LEN)', FFT_LEN));
        Xdb(i,:) = 20*log10(abs(X))';   
    end

    subplot(2,1,1)
    plot(t, real(sig), 'b-', t, imag(sig), 'r-')
    xlim([t(1) t(end)])

    t_range = 1:(length(seg_idx)-1);
    f_range = (((-FFT_LEN/2):(FFT_LEN/2-1))/FFT_LEN)*Fs;

    subplot(2,1,2)
    imagesc(f_range, t_range, Xdb)
    view([270 90])
    colormap('jet') 

end

