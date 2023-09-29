%COMPUTES BINARY SENSOR FIELD

% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);


% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;

% create initial pressure distribution using makeDisc
disc_magnitude = 2; % [Pa]
disc_x_pos = 50;    % [grid points]
disc_y_pos = 20;    % [grid points]
disc_radius = 2;    % [grid points]
disc_1 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 3; % [Pa]
disc_x_pos = 80;    % [grid points]
disc_y_pos = 60;    % [grid points]
disc_radius = 1.5;    % [grid points]
disc_2 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 8; % [Pa]
disc_x_pos = 20;    % [grid points]
disc_y_pos = 35;    % [grid points]
disc_radius = 3;    % [grid points]
disc_3 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

%f(x,y)
source.p0 = disc_1 + disc_2 + disc_3;

% define a CARTESIAN circular sensor
sensor_radius = 4e-3;   % [m]
num_sensor_points = 50;
binary.sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

% define a BINARY sensor mask
sensor_x_pos = Nx/2;        % [grid points]
sensor_y_pos = Ny/2;        % [grid points]
sensor_radius = Nx/2 - 22;  % [grid points]
sensor_arc_angle = 3*pi/2;  % [90 degree sensor gap]
sensor.mask = makeCircle(Nx, Ny, sensor_x_pos, sensor_y_pos, sensor_radius, sensor_arc_angle);

% run the simulation
%rows = number of sensors, columns the data read from those sensors over
%time
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);
% reorder the simulation data < (ONLY WORKS ON 2D DATA - fixes the gaps)
sensor_data_reordered = reorderSensorData(kgrid, sensor, sensor_data);

% plot the simulated sensor data
figure;
imagesc(sensor.mask + source.p0, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;