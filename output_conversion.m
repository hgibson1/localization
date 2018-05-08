Np = 13; % Number of light patterns
Ns = 12; % Number of sensors
LED_states_empty = zeros(Ns, Np);
sensor_readings_empty = zeros(Ns, Np);
LED_states_occupied = zeros(Ns, Np);
sensor_readings_occupied = zeros(Ns, Np);

%% Convert output file into matrix
fid = fopen('Data/5_1/May_1_empty2_154.txt');  % Output readings for empty room
fid2 = fopen('Data/5_1/May_1_Hannah2_154.txt'); % Output readings for occupied room

tline = fgetl(fid); % Read a line from the file as a string
while ischar(tline)
    if ~isempty(tline)
        %disp(tline);
        strC = strsplit(tline, ' '); % Split the string by space
        if length(strC) == 5 && strcmp(char(strC(1)), 'State')
            state_id = str2double(char(strC(3)));
            state_id2 = mod(state_id - 1, Np);
        else
            sensor_id = str2double(char(strC(end - 12)));
            LED_states_empty(sensor_id + 1, state_id + 1) = str2double(char(strC(end - 8)));
            sensor_readings_empty(sensor_id + 1, state_id2 + 1) = str2double(char(strC(end)));
        end
    end
    tline = fgetl(fid); % Read the next line
end
fclose(fid);

tline = fgetl(fid2); % Read a line from the file as a string
while ischar(tline)
    if ~isempty(tline)
        %disp(tline);
        strC = strsplit(tline, ' '); % Split the string by space
        if length(strC) == 5 && strcmp(char(strC(1)), 'State')
            state_id = str2double(char(strC(3)));
            state_id2 = mod(state_id - 1, Np);
        else
            sensor_id = str2double(char(strC(end - 12)));
            LED_states_occupied(sensor_id + 1, state_id + 1) = str2double(char(strC(end - 8)));
            sensor_readings_occupied(sensor_id + 1, state_id2 + 1) = str2double(char(strC(end)));
        end
    end
    tline = fgetl(fid2); % Read the next line
end
fclose(fid2);

%% Calculate light transport matrix
% Empty room
f0 = LED_states_empty(:, 1); % Base light
y0 = sensor_readings_empty(:, 1); % Sensor readings under base light
delta_f_empty = LED_states_empty(:, 2:end) - f0 * ones(1, Np - 1);
delta_y_empty = sensor_readings_empty(:, 2:end) - y0 * ones(1, Np - 1);
A0 = (delta_y_empty * delta_f_empty') / (delta_f_empty * delta_f_empty');

% Occupied room
f0 = LED_states_occupied(:, 1); % Base light
y0 = sensor_readings_occupied(:, 1); % Sensor readings under base light
delta_f = LED_states_occupied(:, 2:end) - f0 * ones(1, Np - 1);
delta_y = sensor_readings_occupied(:, 2:end) - y0 * ones(1, Np - 1);
A = (delta_y * delta_f') / (delta_f * delta_f');
E = A - A0;

save('Data/5_1/Hannah2.mat', 'A0', 'A', 'E');