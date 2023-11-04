clear all output
% This is a 360 Sensor Case Using a Cartesian Circle
% The code performs a forward simulation as well as reconstruction
% This code is pulled **exactly** from k-wave sources

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
disc_x_pos = 64;    % [grid points]
disc_y_pos = 64;    % [grid points]
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

%plotting
sensor_radius_grid_points = round(sensor_radius / kgrid.dx);
binary_sensor_mask = makeCircle(kgrid.Nx, kgrid.Ny, kgrid.Nx/2 + 1, kgrid.Ny/2 + 1, sensor_radius_grid_points, 360);


hold on
% plot the simulated sensor data
% IF YOU WANT TO SEE THE SENSOR MASK, MAKE EQUIVALENT BINARY SENSOR
% MASK AND ADD TO source.p0
figure;
imagesc(source.p0 + binary_sensor_mask, [-1, 1]);
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
img = imagesc(p0_estimate + binary_sensor_mask, [-1, 1]);
maskedRgbImage = img .* cast(~binary_sensor_mask, 'like', img);
maskedRgbImage;
figure;
colormap(getColorMap);
ylabel('x');
xlabel('y');
colorbar;
hold off

clear all output

%% Notes: It seems as if the reconstructed plot is accurate within sensor 
%% boundaries, but outside not so much