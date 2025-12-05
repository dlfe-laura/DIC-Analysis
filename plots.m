%%

%--------------------------------------------------------------------------------------
% Plot the first image + grid boundaries
%-------------------------------------------------------------------------------------

% 1. Display the first image in the series.
% The '[]' is used to scale the intensity range automatically for 'double' type images.
figure;
imshow(ImageSeries{1}, []);

% 2. Hold the current axes to plot additional graphics on top.
hold on;

% 3. Calculate the dimensions required for the rectangle function: [x, y, width, height]
% The (x, y) coordinates for image display usually correspond to the top-left corner.
rect_x = x_min;
rect_y = y_min;
rect_width = x_max - x_min;
rect_height = y_max - y_min;

% 4. Plot the rectangle with the transparent fill
rectangle('Position', [rect_x, rect_y, rect_width, rect_height], ...
          'EdgeColor', 'r', ...             % Boundary color (Red)
          'LineWidth', 2, ...               % Line width
          'FaceColor', [1 0 0], ...         % Fill color (Blue)
          'FaceAlpha', 0.2);                % Transparency level (0.1 for almost transparent)

% 5. Release the hold
hold off;
