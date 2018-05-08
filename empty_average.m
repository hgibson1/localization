Np = 13; % Number of light patterns
Ns = 12; % Number of sensors
LED_states_empty = zeros(Ns, Np);
sensor_readings_empty = zeros(Ns, Np);
LED_states_occupied = zeros(Ns, Np);
sensor_readings_occupied = zeros(Ns, Np);

num_recordings = 6;
A0s = zeros(Ns, Np - 1, num_recordings);

%% Convert output file into matrix
for i = 1:num_recordings
    fid = fopen(['Data/5_3/empty/', num2str(i), 'out.txt']);  % Output readings for empty room
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
    f0 = LED_states_empty(:, 1); % Base light
    y0 = sensor_readings_empty(:, 1); % Sensor readings under base light
    delta_f_empty = LED_states_empty(:, 2:end) - f0 * ones(1, Np - 1);
    delta_y_empty = sensor_readings_empty(:, 2:end) - y0 * ones(1, Np - 1);
    A0s(:, :, i) = (delta_y_empty * delta_f_empty') / (delta_f_empty * delta_f_empty');
end

figure; hold on;
for i = 1:num_recordings
    A0 = A0s(:, :, i);
    plot(A0(:));
end

A0 = mean(A0s, 3);
A = A0;
As = A0s;
save('./Data/5_3/empty.mat', 'As', 'A');