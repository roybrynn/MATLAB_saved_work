% BASIC ITERATIVE SOLUTION - K-WAVE

% THE FOLLOWING CODE IS FROM BASIC_RECON

% 360 Sensor Case
% HAVE TO CLEAR ALL OUTPUT BEFORE USING

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
disc_magnitude = 5; % [Pa]
disc_x_pos = 50;    % [grid points]
disc_y_pos = 50;    % [grid points]
disc_radius = 8;    % [grid points]
disc_1 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 3; % [Pa]
disc_x_pos = 80;    % [grid points]
disc_y_pos = 60;    % [grid points]
disc_radius = 5;    % [grid points]
disc_2 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

source.p0 = disc_1 + disc_2;

% define a centered circular sensor
sensor_radius = 4e-3;   % [m]
num_sensor_points = 50;
cart_sensor_mask = makeCartCircle(sensor_radius, num_sensor_points);
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

hold on
% plot the simulated sensor data
% IF YOU WANT TO SEE THE SENSOR MASK, MAKE EQUIVALENT BINARY SENSOR
% MASK AND ADD TO source.p0
figure;
imagesc(source.p0, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;
hold off

%reset the pressure
source.p0 = 0;

% create the grid, same as before
kgrid_recon = kgrid;

% create a binary sensor mask of an equivalent continuous circle
sensor_radius_grid_points = round(sensor_radius / kgrid_recon.dx);
binary_sensor_mask = makeCircle(kgrid_recon.Nx, kgrid_recon.Ny, kgrid_recon.Nx/2 + 1, kgrid_recon.Ny/2 + 1, sensor_radius_grid_points, 360);

% assign to sensor structure
sensor.mask = binary_sensor_mask;

% interpolate data to remove the gaps and assign to sensor structure
sensor.time_reversal_boundary_data = interpCartData(kgrid_recon, sensor_data, cart_sensor_mask, binary_sensor_mask);

% run the time reversal reconstruction
p0_estimate = kspaceFirstOrder2D(kgrid_recon, medium, source, sensor);

% plot the simulated RECON sensor data
hold on
figure;
imagesc(p0_estimate, [-1, 1]);
colormap(getColorMap);
ylabel('x');
xlabel('y');
colorbar;
hold off

% ALL NEW CODE STARTS BELOW - Iterative Work

% load an image for the initial pressure distribution
p0_image = loadImage('basic_recon_sensor_picture.png');

% make it binary
p0_image = double(p0_image > 0);

% smooth the initial pressure distribution
p0 = smooth(p0_image, true);

% assign to the source structure
source.p0 = p0;

% use the sensor points as sources in time reversal
source.p_mask = sensor.mask;

% time reverse and assign the data
source.p = fliplr(sensor_data);	

% enforce, rather than add, the time-reversed pressure values
source.p_mode = 'dirichlet';    

% set the simulation to record the final image (at t = 0)
sensor.record = {'p_final'};

% run the time reversal reconstruction - did this earlier
%p0_estimate = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% apply a positivity condition
p0_estimate.p_final = p0_estimate.p_final .* (p0_estimate.p_final > 0);

% set the initial pressure to be the latest estimate of p0
source.p0 = p0_estimate.p_final;

% set the simulation to record the time series
sensor = rmfield(sensor, 'record');

% calculate the time series using the latest estimate of p0
sensor_data2 = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% calculate the error in the estimated time series
data_difference = sensor_data - sensor_data2;

% assign the data_difference as a time-reversal source
source.p_mask = sensor.mask;
source.p = fliplr(data_difference);
source = rmfield(source,'p0');
source.p_mode = 'dirichlet';   

% set the simulation to record the final image (at t = 0)
sensor.record = {'p_final'};

% run the time reversal reconstruction
p0_update = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% add the update to the latest image
p0_estimate.p_final = p0_estimate.p_final + p0_update.p_final;

% apply a positivity condition
p0_estimate.p_final = p0_estimate.p_final .* (p0_estimate.p_final > 0);

% plot the MULTIPLE X simulated RECON sensor data
hold on
figure;
imagesc(p0_estimate.p_final + sensor.mask, [-1, 1]);
colormap(getColorMap);
ylabel('x');
xlabel('y');
colorbar;
hold off
clear all output;