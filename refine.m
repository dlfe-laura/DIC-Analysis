%%

%--------------------------------------------------------
% Import images and resize them
% -------------------------------------------------------

%-------------------------------------------------------------------------------
myFolder = 'C:\Users\la3314de\Downloads\drive-download-20251113T090234Z-1-001';
scaleFactor = 0.25;

% Create a list of all the pictures

filePattern = fullfile(myFolder, '*.tif');
imageList = dir(filePattern);

numFiles = length(imageList);

% Pre-allocate a cell array to store the resized images.
ImageSeries = cell(1, numFiles);

% Loop through all found images, resize, and store them in memory
for k = 1:numFiles
    baseFileName = imageList(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
        
    try
        % 1. Read image 
        img = imread(fullFileName);
            
        % 2. Convert to double 
        img = double(img);
            
        % 3. Resize using the specified scale factor and interpolation
        % Target size is 1644 (H) x 1096 (W) based on original 6576x4384 and 0.25 factor.
         resizedImg = imresize(img, scaleFactor, "bicubic");
            
        % 4. Store the resized image in the cell array
        ImageSeries{k} = resizedImg; % This is the list of all RESIZED IMAGES
            
        catch ME
            % Only print if an error occurs
            fprintf('Error processing %s: %s\n', baseFileName, ME.message);
            % Skip this image if it fails
            ImageSeries{k} = []; 
    end
end

% Check resized size
[H, W] = size(ImageSeries{1});

fprintf('Resized image dimensions: [%d × %d]\n', H, W);
%-------------------------------------------------------------------------------

%%

%-------------------------------------------------------
% Define parameters
%--------------------------------------------------------

search_range_px = 2;

ref_x = 10;
ref_y = 10;

x_min = 150;
x_max = 950;

y_min = 550;
y_max = 1200;

node_step = 10;

%%

%-----------------------------------------------------------------
% Other necessary steps before loop correlation
%-----------------------------------------------------------------

% We define the vectors for each direction
x_vect = x_min:node_step:x_max; % = start:step:end
y_vect = y_min:node_step:y_max;

[X, Y] = meshgrid(x_vect, y_vect); % The 2D matrix of the vectors

% The number of nodes in our grid
NUM_nodes_x = length(x_vect);
NUM_nodes_y = length(y_vect);
fprintf('The number of total nodes is %d\n', NUM_nodes_x * NUM_nodes_y);

% We create zero matrices to store the displacement for each node
DX = zeros(NUM_nodes_y, NUM_nodes_x);
DY = zeros(NUM_nodes_y, NUM_nodes_x);

num_frames = numFiles; % Count number of frames
num_steps = num_frames - 2; % Calculate the number of analysis steps

% Pre-allocate cell arrays for history
DX_all = cell(numFiles, 1);
DY_all = cell(numFiles, 1);

% 3. NEW: Set Frame 1 to zero (since it is the reference)
DX_all{1} = DX; 
DY_all{1} = DY;

%%
%-----------------------------------------------------------------
% Reference image
%-----------------------------------------------------------------

img1_ref = ImageSeries{1}; % Therefore all the sequence images will compare their nodes to first image

%%

%-------------------------------------------------------------------
% Loop for image correlation
%-------------------------------------------------------------------

tic

for k = 2:num_steps % Start at 1 until num_steps. It just loops through all the images
    %img1_ref = ImageSeries{k}; % The first image is the reference
    img2_def = ImageSeries{k};  % The second image is the deformed

    %incremental displacement,frame k-1 -> k
    dx_inc = zeros(NUM_nodes_y, NUM_nodes_x);     % incremental dx
    dy_inc = zeros(NUM_nodes_y, NUM_nodes_x);     % incremental dy
    CC_store = zeros(NUM_nodes_y, NUM_nodes_x);   % best CC at each node

    % Loop through the nodes
    for i= 1:NUM_nodes_x % Loop through all the nodes between the x boundaries
        for j= 1:NUM_nodes_y % Loop through all the nodes between the y boundaries

           % Current node
            x_0 = x_vect(i);
            y_0 = y_vect(j);

            % Reference image is centered at (x_0, y_0)
            x_min_0 = x_0 - ref_x;
            x_max_0 = x_0 + ref_x;

            y_min_0 = y_0 - ref_y;
            y_max_0 = y_0 + ref_y;


            zone_0 = img1_ref(y_min_0:y_max_0, x_min_0:x_max_0); % Zone where the reference image is centered (in 2D)

            % Initialize CC_max
            CC_max = -Inf; % Start max correlation for this node at minimum value

            % Reset displacements in the deformed image
            u_change = 0;
            v_change = 0;

            % Loop through the possible displacements
            for u = -search_range_px:search_range_px
                for v = -search_range_px:search_range_px

                    % Boundaries of the displaced image
                    x_min_d = x_min_0 + u;
                    x_max_d = x_max_0 + u;

                    y_min_d = y_min_0 + v;
                    y_max_d = y_max_0 + v;

                    if x_min_d < 1 || x_max_d > W || y_min_d < 1 || y_max_d > H
                        continue;
                    end

                    zone_d = img2_def(y_min_d:y_max_d, x_min_d:x_max_d);

                    % Now we check the CC with the CC_max using NCC

                    % Numerator
                    numerator = sum(sum(zone_0 .* zone_d)); % Matrix inner product of the two images

                    % Denominator
                    denominator = sqrt (sum (sum (zone_0 .^2)) * sum (sum(zone_d .^ 2)));

                    % CC
                    CC = numerator/denominator;

                    % Check for maximum correlation
                    if CC > CC_max
                        CC_max = CC;
                        u_change = u;
                        v_change = v;
                    end

                end % v loop 
            end % u loop

            dx_inc(j,i) = u_change;
            dy_inc(j,i) = v_change;
            CC_store(j, i) = CC_max;

        end % loop over j
    end % loop over i

    % Add the small step to the total displacement
    DX = DX + dx_inc;
    DY = DY + dy_inc;
    
    % Store the new Total in the cell array for Frame k+1
    DX_all{k} = DX;
    DY_all{k} = DY;

    img1_ref = img2_def;

end % Image sequence loop
    
toc

%%

%-----------------------------------------------------------------
% Calculate and plot Displacement
%-----------------------------------------------------------------

Mag_all = cell(numFiles, 1); % create a cell list
mag_min = Inf; % initial boundaries for the max and min magnitudes
mag_max = -Inf;

% Loop through all the files and get the displacement for each frame
for k = 1:numFiles
    % Get the X and Y displacement for this frame
    dx = DX_all{k};
    dy = DY_all{k};
    
    if ~isempty(dx) && ~isempty(dy)
        % Calculate Magnitude
        mag = sqrt(dx.^2 + dy.^2);
        
        % Store it
        Mag_all{k} = mag;
        
        % Update Global Limits
        mag_min = min(mag_min, min(mag(:)));
        mag_max = max(mag_max, max(mag(:)));
    end
end

% Generate the animation 
figure('Name', 'Displacement Magnitude', 'Color', 'w');

fprintf('Animating Magnitude Field...\n');

for k = 1:numFiles
    if isempty(Mag_all{k})
        continue;
    end
    
    % Plot using pcolor
    h = pcolor(X, Y, Mag_all{k});
    
    % Style the plot
    set(h, 'EdgeColor', 'none'); % Remove grid lines for a smooth look
    shading interp;              % Smooths the colors between nodes
    axis equal;                  % Ensures pixels are square
    axis([x_min x_max y_min y_max]); % Lock the view window
    set(gca, 'YDir', 'reverse'); % Flip Y axis to match image coordinates
    
    % Color Map and Limits
    colormap(turbo);               % type of contrast
    colorbar;                    % Show the scale
    if mag_min ~= mag_max
        clim([mag_min mag_max]); % Lock the color scale
    end
    
    % Title
    title(sprintf('Displacement \nStep %d / %d', k, numFiles));
    xlabel('X Position [px]');
    ylabel('Y Position [px]');
    
    % Force update
    drawnow;
    pause(0.5); % Adjust speed of animation here
end

fprintf('Animation Complete.\n');

%%

%-----------------------------------------------------------------
% Calculate and plot Strain
%-----------------------------------------------------------------

%% 
%-----------------------------------------------------------------
% STRAIN FIELD CALCULATION
%-----------------------------------------------------------------
fprintf('Calculating Strain and Global Limits...\n');

% Initialize storage
Exx_all = cell(numFiles,1);
Eyy_all = cell(numFiles,1);
Exy_all = cell(numFiles,1);

% Initialize Global Max/Min for consistent coloring
exx_min_global = inf; exx_max_global = -inf;
eyy_min_global = inf; eyy_max_global = -inf;
exy_min_global = inf; exy_max_global = -inf;
u_max_global = 0;

for k = 1:numFiles
    dxk = DX_all{k};
    dyk = DY_all{k};
    
    if isempty(dxk) || isempty(dyk)
        continue;
    end
    
    % --- 1. Calculate Displacement Magnitude Global Max ---
    u_mag_k = sqrt(dxk.^2 + dyk.^2);
    u_max_global = max(u_max_global, max(u_mag_k(:), [], 'omitnan'));

    % --- 2. Calculate Gradients Manually ---
    % MATLAB gradient returns: [horizontal_derivative, vertical_derivative]
    % We do NOT pass spacing here; we calculate per index (pixel)
    [du_dx_raw, du_dy_raw] = gradient(dxk); 
    [dv_dx_raw, dv_dy_raw] = gradient(dyk); 
    
    % Correct for node_step manually to get physical units (strain)
    du_dx = du_dx_raw / node_step;
    du_dy = du_dy_raw / node_step;
    dv_dx = dv_dx_raw / node_step;
    dv_dy = dv_dy_raw / node_step;

    % --- 3. Calculate Strain Tensor Components ---
    % exx = du/dx
    exx = du_dx;
    % eyy = dv/dy
    eyy = dv_dy;
    % exy = 0.5 * (du/dy + dv/dx)
    exy = 0.5 * (du_dy + dv_dx); 

    % Store data
    Exx_all{k} = exx;
    Eyy_all{k} = eyy;
    Exy_all{k} = exy;

    % --- 4. Update Global Limits for Strains ---
    exx_min_global = min(exx_min_global, min(exx(:), [], 'omitnan'));
    exx_max_global = max(exx_max_global, max(exx(:), [], 'omitnan'));
    
    eyy_min_global = min(eyy_min_global, min(eyy(:), [], 'omitnan'));
    eyy_max_global = max(eyy_max_global, max(eyy(:), [], 'omitnan'));
    
    exy_min_global = min(exy_min_global, min(exy(:), [], 'omitnan'));
    exy_max_global = max(exy_max_global, max(exy(:), [], 'omitnan'));
end

fprintf('Global max |u| = %.3f px\n', u_max_global);
fprintf('Calculations Complete.\n');

%% 
%-----------------------------------------------------------------
% Animation 1: Displacement Magnitude (|u|)
%-----------------------------------------------------------------
figure('Name', 'Displacement Magnitude Animation', 'Color', 'w');

for k = 1:numFiles
    dxk = DX_all{k};
    dyk = DY_all{k};
    
    if isempty(dxk)
        continue;
    end
    
    u_mag_k = sqrt(dxk.^2 + dyk.^2);
    
    pcolor(X, Y, u_mag_k);
    shading interp; % or shading interp for smoothing
    colorbar;
    axis equal; 
    axis([x_min x_max y_min y_max]);
    set(gca, 'YDir', 'reverse'); % Match image coordinates
    
    clim([0 u_max_global]);   % Fixed colour scale based on global max
    colormap(turbo);          % 'turbo' is clearer than 'jet' for magnitude
    
    title(sprintf('|u| (pixels) - Frame %d / %d', k, numFiles));
    xlabel('X Position [px]'); 
    ylabel('Y Position [px]');
    
    drawnow;
    pause(0.1); 
end

%% 
%-----------------------------------------------------------------
% Animation 2: Strain Components (Exx, Eyy, Exy)
%-----------------------------------------------------------------
% Creating one figure with 3 subplots to see all strains simultaneously
figure('Name', 'Strain Fields Animation', 'Color', 'w', 'Position', [50, 50, 1500, 500]);

for k = 1:numFiles
