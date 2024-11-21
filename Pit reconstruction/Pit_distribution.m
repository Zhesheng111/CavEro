clear;
close all;
clc;

%% Step 1: Load and Process Data
filename = 'sphere.mat';  
data = load(filename);    

% Check the structure of the data
if ~isfield(data, 'all_data_points')
    error('The loaded data does not contain the variable all_data_points.');
end

% Display a preview of the data
disp('First 10 rows of data:');
disp(data.all_data_points(1:10, :));  % Show the first 10 rows of data to check for NaN

% Extract coordinates from the data
X = data.all_data_points(:, 1);
Y = data.all_data_points(:, 2);
Z = data.all_data_points(:, 3);

% Check for NaN values in the data
if any(isnan(X)) || any(isnan(Y)) || any(isnan(Z))
    disp('Data contains NaN values, processing NaN values...');
    % Remove data points with NaN
    valid_idx = ~isnan(X) & ~isnan(Y) & ~isnan(Z);
    X = X(valid_idx);
    Y = Y(valid_idx);
    Z = Z(valid_idx);
    
    % Alternatively, interpolate to fill NaN values
    % X = fillmissing(X, 'linear');
    % Y = fillmissing(Y, 'linear');
    % Z = fillmissing(Z, 'linear');
    
    % Display the size of the cleaned data
    disp(['Cleaned data size: ', num2str(length(X)), ' valid points']);
else
    disp('No NaN values in the data');
end

%% Step 2: Plane Fitting and Correction
% Convert X, Y, and Z coordinates to column vectors
X = X(:);
Y = Y(:);
Z = Z(:);

% Construct the design matrix [X, Y, 1] and fit the plane using least squares
A = [X, Y, ones(size(X))];
params = A \ Z;  % Solve the least squares problem

% Extract the plane fitting parameters a, b, c
% a = params(1);  % Coefficient for X
% b = params(2);  % Coefficient for Y
% c = params(3);  % Constant term

a = 0;  % Coefficient for X
b = 0;  % Coefficient for Y
c = 0;  % Constant term

% Display the plane equation
disp(['Fitted plane equation: Z = ', num2str(a), ' * X + ', num2str(b), ' * Y + ', num2str(c)]);

%% Step 3: Grid Data and Interpolation
% Define grid resolution
grid_size = 1000;

% Arrange data into a grid
[xq, yq] = meshgrid(linspace(min(X), max(X), grid_size), ...
                    linspace(min(Y), max(Y), grid_size));

% Interpolate Z data
zq = griddata(X, Y, Z, xq, yq, 'linear');

% Replace NaN values with nearest neighbor interpolation
nan_mask = isnan(zq);
if any(nan_mask(:))
    zq(nan_mask) = griddata(X, Y, Z, xq(nan_mask), yq(nan_mask), 'nearest');
end

%% Step 4: Plot Contour Map with Concentric Circles
figure;
contourf(xq, yq, zq, 20);
colorbar;
hold on;
xlabel('X');
ylabel('Y');
title('Contour Plot with Concentric Circles');

% Set center coordinates and maximum radius
center_x = 879;
center_y = 891;
max_radius = 916.206112751104;
disp(['Maximum radius: ', num2str(max_radius)]);

% Draw concentric circles
num_circles = 20;
radii = linspace(0, max_radius, num_circles);

for i = 1:num_circles
    % Current radius
    r = radii(i);
    
    % Draw concentric circles
    theta = linspace(0, 2 * pi, 100);
    x_circle = center_x + r * cos(theta);
    y_circle = center_y + r * sin(theta);
    plot(x_circle, y_circle, 'w--', 'LineWidth', 1); % White dashed lines for circles
end
hold off;

%% Step 5: Plot 3D Surface and Fitted Plane
figure;
% Plot interpolated 3D surface
surf(xq, yq, zq, 'EdgeColor', 'none');
hold on;

% Plot the fitted plane
Z_plane_grid = a * xq + b * yq + c;  % Compute Z values using the fitted plane equation
surf(xq, yq, Z_plane_grid, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Surface Plot with Fitted Plane');
colorbar;
view(3);
hold off;

%% Step 6: Calculate Volume Below Plane for Each Concentric Ring
% Initialize results for volumes
ring_volumes = zeros(1, num_circles);
volume_per_area = zeros(1, num_circles);       % Volume-to-area ratio for each ring

% Compute grid cell area
dx = xq(1,2) - xq(1,1);
dy = yq(2,1) - yq(1,1);
grid_area = dx * dy;
disp(['Grid cell area: ', num2str(grid_area)]);

for i = 1:num_circles
    % Outer radius of the current ring
    r_outer = radii(i);
    
    if i == 1
        r_inner = 0;
    else
        r_inner = radii(i - 1);
    end
    
    % Create masks for outer and inner circles
    distance_sq = (xq - center_x).^2 + (yq - center_y).^2;
    outer_mask = distance_sq <= r_outer^2;
    inner_mask = distance_sq <= r_inner^2;
    
    % Mask for the current ring = outer circle mask - inner circle mask
    ring_mask = outer_mask & ~inner_mask;
    
    % Get Z values for the current ring region
    Z_plane_current = Z_plane_grid(ring_mask);
    Z_actual = zq(ring_mask);
    
    % Compute height difference (fitted plane - actual Z values)
    z_diff = Z_plane_current - Z_actual;
    
    % Consider only points below the fitted plane (pits)
    z_diff(isnan(z_diff)) = 0;  % Set NaN values to 0
    z_diff(z_diff < 0) = 0;      % Ignore values above the plane
    
    % Compute volume
    ring_volume = sum(z_diff(:)) * grid_area;  % Volume = height difference * cell area
    ring_volumes(i) = ring_volume;
    
    % Compute ring area
    ring_area = pi * (r_outer^2 - r_inner^2);
    
    % Compute volume-to-area ratio
    volume_per_area(i) = ring_volume / ring_area;
    
    % Debug information
    fprintf('Ring %d: Radius = %.2f, Volume = %.4f, Area = %.4f, Volume/Area = %.4f\n', ...
        i, r_outer, ring_volumes(i), ring_area, volume_per_area(i));
end

% Compute total volume
total_volume = sum(ring_volumes);
disp(['Total volume: ', num2str(total_volume)]);

%% Step 7: Output Results and Plot Bar Charts
disp('Volume below the plane for each concentric ring:');
for i = 1:num_circles
    fprintf('Ring %d: Radius = %.2f, Volume = %.4f, Volume/Area = %.4f\n', ...
        i, radii(i), ring_volumes(i), volume_per_area(i));
end
fprintf('Total volume: %.4f\n', total_volume);

% Plot volume-to-area ratio bar chart
figure;
bar(radii, volume_per_area, 'FaceAlpha', 0.7, 'EdgeColor', 'k');
xlabel('Radius of Concentric Circles');
ylabel('Unit Volume (Volume per Area)');
title('Unit Volume below Fitted Plane within Concentric Rings');
grid on;

% Add total volume-to-area ratio label
total_area = pi * max_radius^2;
text(max(radii)*0.6, max(volume_per_area)*0.9, ...
     ['Total Volume/Total Area = ', num2str(total_volume / total_area, '%.4f')], ...
     'FontSize', 12, 'Color', 'b', 'BackgroundColor', 'w');

% Plot actual volume bar chart
figure;
bar(radii, ring_volumes, 'FaceAlpha', 0.7, 'EdgeColor', 'k');
xlabel('Radius of Concentric Circles');
ylabel('Volume below Plane');
title('Volume below Fitted Plane for Concentric Rings');
grid on;

%% Step 8: Save Results to MAT File
% Create a matrix with three columns: radius, volume, volume-to-area ratio
output_data = [radii(:), ring_volumes(:), volume_per_area(:)];

% Save data to MAT file
save('volume_per_ring_data.mat', 'output_data');

% Notify user that the data has been saved
disp('Volume data has been saved to volume_per_ring_data.mat.');
