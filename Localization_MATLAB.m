%clear all; close all; 
clc;

%% Parameters
Lx = 2.800;  % Room length in x dimension (Unit: meter or ft)
Ly = 2.000;  % Room length in y dimension (Unit: meter or ft)
H = 2.710;  % Room height (Unit: meter or ft)

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

%% Localization algorithm
stepsize = Lx / 50;
fprintf(1, 'Generating matrix C...\n');
if exist('C', 'var') ~= 1
    [C, WW, LL] = C_gen_RPI(Lx, Ly, H, 1e-1, stepsize, LED_pos, sensor_pos, -2);
end
fprintf(1, 'Running localization algorithm...\n');
place_ind = 3;
load Data/5_2/101ms/empty.mat;
load(['Data/5_2/101ms/Timmy_place', num2str(place_ind), '.mat']);
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
for k = 1:size(Rec, 1)
    if Rec(k,5) > 0
        rectangle('Position', Rec(k, 1:4), 'FaceColor', [1, 1-Rec(k, 5), 1-Rec(k, 5)], 'LineStyle', 'none');
    else
        rectangle('Position', Rec(k, 1:4), 'FaceColor', [1+Rec(k, 5), 1, 1+Rec(k, 5)], 'LineStyle', 'none');
    end
end
% Plot LEDs
for j = 1:Nf
    h3 = plot(LED_pos(1, j), LED_pos(2, j), 'ko', 'MarkerSize', 8);
end
% Plot sensors
for j = 1:Ns
    h4 = plot(sensor_pos(1, j), sensor_pos(2, j), 'ks', 'MarkerSize', 8);
end
% Plot estimated and ground truth locations
h1 = plot(loc_est(1), loc_est(2), 'kx', 'MarkerSize', 8);
h2 = plot(gt(1), gt(2), 'k.', 'MarkerSize', 20);
loc_error = norm(loc_est - gt, 2);
legend([h1, h2, h3, h4], 'Estimated location', 'Ground truth', 'LEDs', 'sensors');
% Plot the floor bound
rectangle('Position', [0, 0, Lx, Ly], 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1);
set(gca, 'XLim', [0, Lx]); set(gca, 'YLim', [0, Ly]);
xlabel('x (m)'); ylabel('y (m)'); 
title(sprintf('Localization error = %.3fm', loc_error));
fprintf(1, 'Estimated location: (%.4f, %.4f)\n', loc_est(1), loc_est(2));
fprintf(1, 'Localization error = %.3fm\n', loc_error);
