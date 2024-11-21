clear; clc; close all;

img = imread('0.5-3.tif'); 

% Get the size of the image
[img_height, img_width] = size(img);

% Assume a scale factor of 1 pixel = 0.5 microns
scale_factor = 1920/1200;  % Actual physical size per pixel, in microns

% Convert the image to binary to separate the pit regions (assuming white is the pit)
bw_img = img > 128;  % Adjust threshold based on the image

% Fill closed pits
filled_pits = imfill(bw_img, 'holes');

% Remove small noise areas, keep larger pits
min_area_pixels = 5;  % Set a reasonable minimum area threshold, adjust based on the image
filtered_pits = bwareaopen(filled_pits, min_area_pixels);

% Get the boundaries of each pit (using boundary tracing)
boundaries = bwboundaries(filtered_pits);

% Extract the edge image and save it
edge_img = false(size(img));  % Initialize edge image

% Traverse each boundary and add it to the edge image
for k = 1:length(boundaries)
    boundary = boundaries{k};
    edge_mask = poly2mask(boundary(:,2), boundary(:,1), size(img, 1), size(img, 2));
    edge_img = edge_img | edge_mask;  % Use logical OR to combine all edges
end

% Save the extracted edge image
imwrite(edge_img, 'extracted_edges.png');  

% -----------------------------
% Image scale conversion (pixel to actual coordinates)
% Convert image width and height (in pixels) to real-world coordinate system
[X_plane, Y_plane] = meshgrid(1:img_width, 1:img_height);  % Use image pixel size
Z_plane = zeros(size(X_plane));  % Z plane at z=0

% Convert pixel coordinates to real-world coordinates (in microns)
X_plane = X_plane * scale_factor;
Y_plane = Y_plane * scale_factor;

% Add the z=0 plane data (X_plane, Y_plane, Z_plane) to all data points
plane_data_points = [X_plane(:), Y_plane(:), Z_plane(:)];

% ---------------------------
% Traverse each pit boundary, calculate equivalent circle and process ellipsoidal surface data

all_data_points = [];  % Initialize data points array

for k = 1:length(boundaries)
    boundary = boundaries{k};
    
    % Create edge image to extract the area of the edge region
    edge_mask = poly2mask(boundary(:,2), boundary(:,1), size(img, 1), size(img, 2));
    
    % Calculate the area enclosed by the edge (in pixels)
    area_pixels = sum(edge_mask(:));
    
    % Calculate the radius of the equivalent circle, in pixels
    radius_pixels = sqrt(area_pixels / pi);
    
    % Calculate the centroid of the edge
    stats = regionprops(edge_mask, 'Centroid');
    if ~isempty(stats) && isfield(stats, 'Centroid')
        centroid = stats.Centroid;
    else
        warning('No valid region or centroid information detected, skipping this pit boundary');
        continue;  % Skip current loop and process the next pit
    end

    % ---------------------------
    % Generate ellipsoid (adjust Z values based on the given formula)
    
    % Define ellipsoid grid parameters
    [theta, phi] = meshgrid(linspace(0, 2*pi, 100), linspace(0, pi, 100));
    
    % Ellipsoid major and minor axes (defined as 1 here, adjust as needed)
    a = radius_pixels * scale_factor;  % Major axis
    b = radius_pixels * scale_factor;  % Minor axis
    c = @(r) -0.0004 * r^2 + 0.1244 * r + 0.0310;  % Depth function of the ellipsoid
    
    % Calculate ellipsoid coordinates
    X = a * sin(phi) .* cos(theta) + centroid(1) * scale_factor;
    Y = b * sin(phi) .* sin(theta) + centroid(2) * scale_factor;
    
    % Calculate Z coordinates, apply depth function, adjust Z values
    Z = c(radius_pixels) * cos(phi);  % Adjust Z coordinates using the depth relation
    
    % Keep only the part where z <= 0 of the ellipsoid
    Z(Z > 0) = NaN;  % Set z > 0 to NaN, removing the upper half
    
    % Get data of the lower part of the ellipsoid
    ellipsoid_data = [X(:), Y(:), Z(:)];
    all_data_points = [all_data_points; ellipsoid_data];  % Add ellipsoid data points to the all_data_points array
end

% ---------------------------

set(gcf, 'Renderer', 'painters');  

% Plot the horizontal plane and ellipsoid surface (plot them all at once)
figure;
hold on;
axis equal;
view(3);  % Set 3D view
xlabel('X (micrometers)');
ylabel('Y (micrometers)');
zlabel('Z (micrometers)');

% ---------------------------
% Plot the z=0 plane
surf(X_plane, Y_plane, Z_plane, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plot all ellipsoid surfaces
for k = 1:length(boundaries)
    boundary = boundaries{k};
    
    % Create edge image to extract the area of the edge region
    edge_mask = poly2mask(boundary(:,2), boundary(:,1), size(img, 1), size(img, 2));
    
    % Calculate the area enclosed by the edge (in pixels)
    area_pixels = sum(edge_mask(:));
    
    % Calculate the radius of the equivalent circle, in pixels
    radius_pixels = sqrt(area_pixels / pi);
    
    % Calculate the centroid of the edge
    stats = regionprops(edge_mask, 'Centroid');
    if ~isempty(stats) && isfield(stats, 'Centroid')
        centroid = stats.Centroid;
    else
        warning('No valid region or centroid information detected, skipping this pit boundary');
        continue;  % Skip current loop and process the next pit
    end

    % ---------------------------
    % Plot ellipsoid surface
    [theta, phi] = meshgrid(linspace(0, 2*pi, 100), linspace(0, pi, 100));
    
    % Ellipsoid major and minor axes (defined as 1 here, adjust as needed)
    a = radius_pixels * scale_factor;  % Major axis
    b = radius_pixels * scale_factor;  % Minor axis
    c = @(r) -0.0004 * r^2 + 0.1316 * r + 0.0121;  % Depth function of the ellipsoid
    
    % Calculate ellipsoid coordinates
    X = a * sin(phi) .* cos(theta) + centroid(1) * scale_factor;
    Y = b * sin(phi) .* sin(theta) + centroid(2) * scale_factor;
    
    % Calculate Z coordinates, apply depth function, adjust Z values
    Z = c(radius_pixels) * cos(phi);  
    
    % Keep only the part where z <= 0 of the ellipsoid
    Z(Z > 0) = NaN; 
    
    % Plot the ellipsoid surface
    surf(X, Y, Z, 'FaceColor', 'g', 'EdgeColor', 'none');
end

hold off;
title('3D Image: Equivalent Circles and Ellipsoid Surfaces');

% -----------------------------
% Show 3D data points: Display all points below the ellipsoid surface
figure;
scatter3(all_data_points(:, 1), all_data_points(:, 2), all_data_points(:, 3), ...
    2, 'filled');
title('3D Data Points: Points Below All Ellipsoid Surfaces');
xlabel('X (micrometers)');
ylabel('Y (micrometers)');
zlabel('Z (micrometers)');

% -----------------------------
% Merge data: Ellipsoid data + plane data
all_data_points_with_plane = [all_data_points; plane_data_points];

% Export 3D data (ellipsoid surface + horizontal plane) in real-world coordinates (microns)
save('ellipsoid_with_plane_data.mat', 'all_data_points_with_plane');
