clear; clc; close all;

% Read the image
img = imread('0.5-3.tif'); 

% Convert the image to a binary image to isolate pit regions (assuming white is the pit)
bw_img = img > 128;  

% Identify closed pit regions
filled_pits = imfill(bw_img, 'holes');

% Remove small noise regions, keeping only larger pit areas
min_area_pixels = 5;  
filtered_pits = bwareaopen(filled_pits, min_area_pixels);

% Get the edges of each pit (using boundary tracing)
boundaries = bwboundaries(filtered_pits);

% Extract edges
edge_img_no_fill = false(size(img)); 

% Iterate through each boundary and add it to the edge image
for k = 1:length(boundaries)
    boundary = boundaries{k};
    for i = 1:size(boundary, 1)
        edge_img_no_fill(boundary(i, 1), boundary(i, 2)) = 1; 
    end
end

% Save the extracted edge image without filling
imwrite(edge_img_no_fill, 'extracted_edges_no_fill.png');  

% Display the extracted edge image without filling
figure;
imshow(edge_img_no_fill);
title('Pit Edges');

% Extract filled edge image
edge_img_fill = false(size(img));  % Initialize the filled edge image

% Combine all filled edges
for k = 1:length(boundaries)
    boundary = boundaries{k};
    edge_mask = poly2mask(boundary(:,2), boundary(:,1), size(img, 1), size(img, 2));
    edge_img_fill = edge_img_fill | edge_mask;  % Combine edges
end

% Save the extracted filled edge image
imwrite(edge_img_fill, 'extracted_edges_fill.png');  

% Display the extracted filled edge image
figure;
imshow(edge_img_fill);
title('Filled Edges');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Plot the image with equivalent circles
figure;
imshow(img); hold on;

% Image scaling factor
scale_factor = 1920/1200;  

% Iterate through each pit boundary, calculate the equivalent circle, and plot it
for k = 1:length(boundaries)
    boundary = boundaries{k};
    
    plot(boundary(:,2), boundary(:,1), 'b-', 'LineWidth', 0.5);  

    % Create an edge mask to extract the region enclosed by the boundary
    edge_mask = poly2mask(boundary(:,2), boundary(:,1), size(img, 1), size(img, 2));
    
    % Calculate the area enclosed by the boundary (in pixels)
    area_pixels = sum(edge_mask(:));
    
    % Compute the radius of the equivalent circle
    radius_pixels = sqrt(area_pixels / pi);
    
    % Calculate the centroid of the boundary
    stats = regionprops(edge_mask, 'Centroid');
    centroid = stats.Centroid;

    % Draw the equivalent circle
    viscircles(centroid, radius_pixels, 'EdgeColor', 'r', 'LineWidth', 0.5);
end

hold off; 
title('Equivalent Circles of Pits');
