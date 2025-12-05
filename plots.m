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

%%
%-------------------------------------------------------
% Save Displacement Vector Animation as GIF
%-------------------------------------------------------

% 1. Define Output Filename
gifFilename = 'VectorField_Animation.gif';

% 2. Setup Figure
hFig = figure('Name', 'Displacement Vector Field Animation');
% Optional: Set specific figure size so the GIF resolution is consistent
set(hFig, 'Position', [100, 100, 800, 600]); 

scale_factor = 1; 

% Flag to track the first valid frame for writing mode
firstFrame = true;

% 3. Loop through time steps (Removed the 'while' loop)
for k = 1:length(DX_all)
    
    % Skip empty frames
    if isempty(DX_all{k}) 
        continue; 
    end
    
    % --- Plotting Routine (Same as your code) ---
    imagesc(ImageSeries{k+1}); 
    colormap gray; 
    axis image; 
    axis off;   
    hold on;
    
    q = quiver(X, Y, ...
               DX_all{k} * scale_factor, ...
               DY_all{k} * scale_factor, ...
               'r', 'LineWidth', 1.5);
    
    q.AutoScale = 'off'; 
    title(sprintf('Frame: %d | Total Displacement', k));
    hold off;
    
    % Force the update so getframe captures the correct state
    drawnow; 
    
    % --- GIF Capture Routine ---
    
    % 1. Capture the current figure as an image
    frame = getframe(hFig); 
    im = frame2im(frame); 
    
    % 2. Convert RGB image to indexed image (required for GIF)
    [imind, cm] = rgb2ind(im, 256); 
    
    % 3. Write to the GIF file
    if firstFrame
        % For the first frame, create the file with LoopCount = inf (infinite loop)
        imwrite(imind, cm, gifFilename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
        firstFrame = false;
    else
        % For subsequent frames, append to the existing file
        imwrite(imind, cm, gifFilename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end
end

fprintf('Animation saved as %s\n', gifFilename);
