%clear all; close all; 
clc;

%% Parameters
Lx = 2.800;  % Room length in x dimension (Unit: meter or ft)
Ly = 2.000;  % Room length in y dimension (Unit: meter or ft)
H = 2.710;  % Room height (Unit: meter or ft)
place_ind = 1; % Index of ground truth location (1 - 3)

% Parameters for the localization algorithm
sigma = 1e-8; % weight of the l_2 penalty term
thr = 0.6;   % threshold

%% Source and sensor locations
% Position of each LED and sensor, each column corresponds to an LED/sensor
% The positions are the x, y, z measurements in a Cartesian coordinate
% Choose a corner of the floor as the origin when measuring
% Please fill in new measurements in meter or ft. 
% Column #1 - #12 should correspond to LED/sensor #1 - #12 respectively.
LED_pos = [0.380, 1.010, 1.650, 2.260, 0.385, 1.018, 1.663, 2.266, 0.373, 1.022, 1.661, 2.262;
           0.390, 0.390, 0.390, 0.390, 1.000, 1.000, 1.000, 1.000, 1.610, 1.610, 1.610, 1.610;
           H, H, H, H, H, H, H, H, H, H, H, H];
sensor_pos = [0.505, 1.135, 1.775, 2.385, 0.510, 1.143, 1.788, 2.391, 0.498, 1.147, 1.786, 2.387;
              0.390, 0.390, 0.390, 0.390, 1.000, 1.000, 1.000, 1.000, 1.610, 1.610, 1.610, 1.610;
              H, H, H, H, H, H, H, H, H, H, H, H];
Nf = size(LED_pos, 2);
Ns = size(sensor_pos, 2);

%% Read output files
LED_states_empty = zeros(Ns, Nf + 1);
sensor_readings_empty = zeros(Ns, Nf + 1);
LED_states_occupied = zeros(Ns, Nf + 1);
sensor_readings_occupied = zeros(Ns, Nf + 1);

list = ls('./Data/empty/*out.txt');
A0s = [];
for i = 1:size(list, 1)
    timestamp = str2double(cell2mat(regexp(list(i, :),'\d*','Match')));
    fid = fopen(['./Data/empty/', num2str(timestamp), 'out.txt']);
    % Check whether file is complete
    states_flag = zeros(1, Nf + 1);
    tline = fgetl(fid);
    while ischar(tline)
        if ~isempty(tline)
            strC = strsplit(tline, ' '); % Split the string by space
            if length(strC) == 5 && strcmp(char(strC(1)), 'State')
                state_id = str2double(char(strC(3)));
                state_id2 = mod(state_id - 1, Nf + 1);
                if state_id >= 0 && state_id < Nf + 1
                    states_flag(state_id + 1) = 1;
                end
            else
                sensor_id = str2double(char(strC(end - 12)));
                if sensor_id >= 0 && sensor_id < Ns && state_id >= 0 && state_id < Nf + 1
                    LED_states_empty(sensor_id + 1, state_id + 1) = str2double(char(strC(end - 8)));
                    sensor_readings_empty(sensor_id + 1, state_id2 + 1) = str2double(char(strC(end)));
                end
            end
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    if length(find(states_flag == 1)) == Nf + 1
        fprintf(1, '%d ', timestamp);
        f0 = LED_states_empty(:, 1); % Base light
        y0 = sensor_readings_empty(:, 1); % Sensor readings under base light
        delta_f = LED_states_empty(:, 2:end) - f0 * ones(1, Nf);
        delta_y = sensor_readings_empty(:, 2:end) - y0 * ones(1, Nf);
        A0 = (delta_y * delta_f') / (delta_f * delta_f');
        A0s = cat(3, A0s, A0);
    end
end
fprintf(1, '\n');
list = ls('./Data/occupied/*out.txt');
As = [];
for i = 1:size(list, 1)
    timestamp = str2double(cell2mat(regexp(list(i, :),'\d*','Match')));
    fid = fopen(['./Data/occupied/', num2str(timestamp), 'out.txt']);
    % Check whether file is complete
    states_flag = zeros(1, Nf + 1);
    tline = fgetl(fid);
    while ischar(tline)
        if ~isempty(tline)
            strC = strsplit(tline, ' '); % Split the string by space
            if length(strC) == 5 && strcmp(char(strC(1)), 'State')
                state_id = str2double(char(strC(3)));
                state_id2 = mod(state_id - 1, Nf + 1);
                if state_id >= 0 && state_id < Nf + 1
                    states_flag(state_id + 1) = 1;
                end
            else
                sensor_id = str2double(char(strC(end - 12)));
                if sensor_id >= 0 && sensor_id < Ns && state_id >= 0 && state_id < Nf + 1
                    LED_states_occupied(sensor_id + 1, state_id + 1) = str2double(char(strC(end - 8)));
                    sensor_readings_occupied(sensor_id + 1, state_id2 + 1) = str2double(char(strC(end)));
                end
            end
        end
        tline = fgetl(fid);
    end
    fclose(fid);
    if length(find(states_flag == 1)) == Nf + 1
        fprintf(1, '%d ', timestamp);
        f0 = LED_states_occupied(:, 1); % Base light
        y0 = sensor_readings_occupied(:, 1); % Sensor readings under base light
        delta_f = LED_states_occupied(:, 2:end) - f0 * ones(1, Nf);
        delta_y = sensor_readings_occupied(:, 2:end) - y0 * ones(1, Nf);
        A = (delta_y * delta_f') / (delta_f * delta_f');
        As = cat(3, As, A);
    end
end
fprintf(1, '\n');
A0 = mean(A0s, 3);
A = mean(As, 3);

%% Localization algorithm
stepsize = Lx / 50;
fprintf(1, 'Generating matrix C...\n');
if exist('C', 'var') ~= 1
    [C, WW, LL] = C_gen_RPI(Lx, Ly, H, 1e-1, stepsize, LED_pos, sensor_pos, -2);
end
fprintf(1, 'Running localization algorithm...\n');

gts = [0.535, 0.52;
       1.52, 1.125;
       2.34, 0.48];
gt = gts(place_ind, :);
E = A - A0;
iterations = 2;
%[delta_alphas, shadows, ~, delta_alpha2, loc_est] = loc_algo_3D_v2(E, C, WW, LL, thr, sigmas, stepsize, 1);
[delta_alpha, locs, ~, ~, ~, ~] = multistep_v2_negOnly(E(:), C, WW, LL, thr, sigma * ones(1, 2), stepsize, iterations); loc_est = locs(:, iterations)';
fprintf(1, 'Complete!\n');

%% Plotting the albedo change map and estimated location
[~, ~, delta_alpha2_2D] = plot_alpha(delta_alpha, WW, LL, stepsize);
Rec = DrawRectangle2(delta_alpha2_2D, WW, LL, stepsize);

figure; hold on; axis equal;
set(gca, 'YDir','reverse');
for k = 1:size(Rec, 1)
    if Rec(k,5) > 0
        rectangle('Position', Rec(k, [2, 1, 3, 4]), 'FaceColor', [1, 1-Rec(k, 5), 1-Rec(k, 5)], 'LineStyle', 'none');
    else
        rectangle('Position', Rec(k, [2, 1, 3, 4]), 'FaceColor', [1+Rec(k, 5), 1, 1+Rec(k, 5)], 'LineStyle', 'none');
    end
end
% Plot LEDs
for j = 1:Nf
    h3 = plot(LED_pos(2, j), LED_pos(1, j), 'ko', 'MarkerSize', 8);
end
% Plot sensors
for j = 1:Ns
    h4 = plot(sensor_pos(2, j), sensor_pos(1, j), 'ks', 'MarkerSize', 8);
end
% Plot estimated and ground truth locations
h1 = plot(loc_est(2), loc_est(1), 'kx', 'MarkerSize', 8);
h2 = plot(gt(2), gt(1), 'k.', 'MarkerSize', 20);
loc_error = norm(loc_est - gt, 2);
legend([h1, h2, h3, h4], 'Estimated location', 'Ground truth', 'LEDs', 'sensors');
% Plot the floor bound
rectangle('Position', [0, 0, Ly, Lx], 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1);
set(gca, 'XLim', [0, Ly]); set(gca, 'YLim', [0, Lx]);
xlabel('y (m)'); ylabel('x (m)'); 
title(sprintf('Localization error = %.3fm', loc_error));
fprintf(1, 'Estimated location: (%.4f, %.4f)\n', loc_est(1), loc_est(2));
fprintf(1, 'Localization error = %.3fm\n', loc_error);
