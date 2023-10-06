%Clear all output

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
%kgrid.makeTime(medium.sound_speed);

% create initial pressure distribution using makeDisc (MAYBE CHANGE THESE
% AROUND)
disc_magnitude = 2; % [Pa]
disc_x_pos = 50;    % [grid points]
disc_y_pos = 70;    % [grid points]
disc_radius = 2;    % [grid points]
disc_1 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 0; % [Pa]
disc_x_pos = 40;    % [grid points]
disc_y_pos = 25;    % [grid points]
disc_radius = 1.5;    % [grid points]
disc_2 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

disc_magnitude = 0; % [Pa]
disc_x_pos = 20;    % [grid points]
disc_y_pos = 35;    % [grid points]
disc_radius = 3;    % [grid points]
disc_3 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

%f(x,y)
source.p0 = disc_1 + disc_2 + disc_3;

% define a centered circular sensor
sensor_radius = 4e-3;   % [m]
num_sensor_points = 50;
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);


%hold on
%figure;
%imagesc(source.p0  + sensor.mask, [-1, 1]);
%imagesc(sensor.mask, [-1, 1]);
%colormap(getColorMap);
%ylabel('Sensor Position');
%xlabel('Time Step');
%colorbar;
%hold off

% THE FOLLOWING CODE COMES FROM:
%http://www.medphys.ucl.ac.uk/research/mle/pdf_files/J_2010_Treeby_JBO_k-Wave_MATLAB_toolbox.pdf

% create the time array
kgrid.t_array = makeTime(kgrid, medium.sound_speed);
% run the forward simulation

sensor_data = kspaceFirstOrder2D(kgrid,medium, source, sensor);

% reorder the simulation data < (ONLY WORKS ON 2D DATA - fixes the gaps)
%sensor_data = reorderSensorData(kgrid, sensor, sensor_data); %this doesnt seem to solve the problem


% reset the initial pressure
source.p0 = 0;
% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data;
% run the time reversal reconstruction
p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% NOW BACK TO KWAVE (kspaceLineRecon)

% define a binary line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;

% plot the simulated RECON sensor data
hold on
figure;
imagesc(p0_recon, [-1, 1]);
colormap(getColorMap);
ylabel('x');
xlabel('y');
colorbar;

% reconstruct the initial pressure
p_xy = kspaceLineRecon(sensor_data.', dy, kgrid.dt, medium.sound_speed, 'Plot', true);
hold off
clear all output