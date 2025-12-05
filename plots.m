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

%%

%-------------------------------------------------
% Displacement Field Heat map
%-----------------------------------------------------

%-------------------------------------------------------
% 1. Pre-calculate Global Max for Color Consistency
%-------------------------------------------------------
max_disp = 0;
for k = 1:length(DX_all)
    if isempty(DX_all{k}), continue; end
    
    curr_mag = sqrt(DX_all{k}.^2 + DY_all{k}.^2); 
    current_max = max(curr_mag(:)); 
    
    if ~isnan(current_max)
        max_disp = max(max_disp, current_max);
    end
end

% Safety check for max_disp
if max_disp <= 1e-6 
    fprintf('Warning: Max displacement is ~0. Setting scale to 1.\n');
    max_disp = 1; 
end

%-------------------------------------------------------
% 2. Setup Figure and GIF Parameters
%-------------------------------------------------------
gifFilename = 'Displacement_Magnitude.gif';

hFig = figure;
set(hFig, 'Name', 'Displacement Only');
% Set a fixed position so the GIF resolution remains constant
set(hFig, 'Position', [100, 100, 800, 600]); 

firstFrame = true; % Flag to handle the first write operation

%-------------------------------------------------------
% 3. Animation and Saving Loop
%-------------------------------------------------------
% Note: We removed the 'while' loop. We only iterate through the data once to save.
for k = 1:length(DX_all)
    
    % Skip empty frames
    if isempty(DX_all{k})
        continue; 
    end
    
    % --- Calculation & Plotting ---
    Mag = sqrt(DX_all{k}.^2 + DY_all{k}.^2); 
    
    pcolor(X, Y, Mag);
    shading interp;
    colormap jet;
    axis equal; axis tight;
    
    % Lock color scale bar to global max
    clim([0 max_disp]); 
    
    hBar = colorbar;
    ylabel(hBar, 'Displacement [px]');
    title(sprintf('Frame %d | Max Disp: %.4f', k, max(Mag(:))));
    
    % Update figure drawing
    drawnow;
    
    % --- GIF Capture Routine ---
    
    % 1. Capture the figure content
    frame = getframe(hFig); 
    im = frame2im(frame); 
    
    % 2. Convert to indexed image (256 colors)
    [imind, cm] = rgb2ind(im, 256); 
    
    % 3. Write to GIF
    if firstFrame
        % Create the file. 'Loopcount', inf makes it repeat forever.
        imwrite(imind, cm, gifFilename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
        firstFrame = false;
    else
        % Append subsequent frames to the file
        imwrite(imind, cm, gifFilename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
    end

    %%

end

fprintf('Animation successfully saved as %s\n', gifFilename);

%%
%----------------------------------------------------------------------
% Strain plots - GIF Generation
%-----------------------------------------------------------------------

% 1. SETUP OUTPUT FILE
gif_filename = 'Strain_Animation.gif'; % Name of the file to save
delay_time = 0.5;                      % Time per frame (seconds)

% Pre-Calculate global limits (Same as your original code)
get_limits = @(data) [min(data(:)) max(data(:))];
clim_x  = get_limits(Epsilon_x);
clim_y  = get_limits(Epsilon_y);
clim_xy = get_limits(Gamma_xy);

% Safety Check
if clim_x(1) == clim_x(2),  clim_x(2)  = clim_x(2) + 1e-6; end
if clim_y(1) == clim_y(2),  clim_y(2)  = clim_y(2) + 1e-6; end
if clim_xy(1) == clim_xy(2), clim_xy(2) = clim_xy(2) + 1e-6; end

% Create Figure
hFig = figure('Units', 'normalized', 'Position', [0.1 0.3 0.8 0.4]); 
set(hFig, 'Name', 'Strain Field Animation', 'Color', 'w'); % Set background to white for cleaner GIF

fprintf('Generating GIF... Please wait.\n');

% 2. ANIMATION LOOP (Run through steps exactly once)
for k = 1:num_steps
    
    % --- Plot 1: Epsilon X ---
    subplot(1, 3, 1);
    pcolor(X, Y, Epsilon_x(:,:,k));
    shading interp; axis equal; axis tight;
    colormap(gca, jet); 
    clim(clim_x);      
    colorbar;
    title(sprintf('E_{xx} (Normal X)\nFrame: %d', k));

    % --- Plot 2: Epsilon Y ---
    subplot(1, 3, 2);
    pcolor(X, Y, Epsilon_y(:,:,k));
    shading interp; axis equal; axis tight;
    colormap(gca, jet); 
    clim(clim_y);      
    colorbar;
    title(sprintf('E_{yy} (Normal Y)\nFrame: %d', k));

    % --- Plot 3: Gamma XY ---
    subplot(1, 3, 3);
    pcolor(X, Y, Gamma_xy(:,:,k));
    shading interp; axis equal; axis tight;
    colormap(gca, jet); 
    clim(clim_xy);     
    colorbar;
    title(sprintf('\\gamma_{xy} (Shear)\nFrame: %d', k));
    
    % Force the update so getframe captures the current state
    drawnow; 
    
    % 3. CAPTURE AND WRITE FRAME
    frame = getframe(hFig); 
    im = frame2im(frame); 
    [imind, cm] = rgb2ind(im, 256); 
    
    % Write to the GIF file
    if k == 1 
        % For the first frame, create the file and set LoopCount to infinity
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time); 
    else 
        % For subsequent frames, append to the existing file
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time); 
    end 
end

fprintf('Done! Saved as %s\n', gif_filename);

