%% 
%-------------------------------------------------------------------------
% First we define the parameters
%-------------------------------------------------------------------------

    % Node points. We define a rectangular grid with limits nd step sizes
    x_min = 41;
    x_max = 1056;
    node_step_x = 20; % the step size between analysis node

    y_min = 41;
    y_max = 1604;
    node_step_y = 20;
    
    % We define the vectors for each direction to define the 2D-matrices
    x_vect = x_min:node_step_x:x_max; % = start:step:end
    y_vect = y_min:node_step_y:y_max;
    
    % The number of nodes
    NUM_nodes_x = length(x_vect);
    NUM_nodes_y = length(y_vect);

    %The 2D matrices
    [X, Y] = meshgrid(x_vect, y_vect);

    % Search window size (smaller window that moves around inside my analysis window to perform a local calculation at each grid node.)
    % It defines the range of possible integer displacements (u, v) around node x_i
    search_range_px = 16;
    
    % Correlation window size (the specific size of a pattern I am looking for). 
    % The template will slide across my entire analysis window to find matches.
    corr_half = 15;
    ref_x= 41 ;
    ref_y = 41 ;
%% 

%-------------------------------------------------------------------------
% Read, plot image, check how many pixels it has, adjust my parameters and
% limits and then resize the images. After start looping
%-------------------------------------------------------------------------

    % --- Setup ---
    myFolder = 'C:\Users\la3314de\Downloads\drive-download-20251113T090234Z-1-001'; 
    scaleFactor = 0.25;
    
    % Get a list of all .tif files in the folder
    filePattern = fullfile(myFolder, '*.tif');
    imageList = dir(filePattern);
    
    % Check if files were found
    if isempty(imageList)
        warning('No TIF files found in this folder: %s', myFolder);
        return;
    end
    
    numFiles = length(imageList);
    
    % Pre-allocate a cell array to store the resized images.
    ImageSeries = cell(1, numFiles);
    
    fprintf('Found %d TIF images. Starting series import and resize...\n', numFiles);

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

    % Verify the first resized image size (Optional check for debugging)
    if ~isempty(ImageSeries{1})
        [H, W] = size(ImageSeries{1});
        fprintf('Import and resize complete. First image size: [%d x %d].\n', H, W);
    else
        % If the first image failed to load, stop
        warning('First image failed to load or resize.');
        return;
    end
%% 
%-------------------------------------------------------------------------
% Define and store results
%-------------------------------------------------------------------------

    num_frames = numFiles; % Count number of frames

    num_steps = num_frames - 1; % Calculate the number of analysis steps

    % Prealocate matrices to store the results
    U_int_disp = zeros(NUM_nodes_y, NUM_nodes_x, num_frames - 1); % Integer displacement in X
    V_int_disp = zeros(NUM_nodes_y, NUM_nodes_x, num_frames - 1); % Integer displacement in Y
    CC_max_store = zeros(NUM_nodes_y, NUM_nodes_x, num_frames - 1); % Max correlation coefficient

%-------------------------------------------------------------------------
% Loop over DIC analysis nodes (grid) 
%-------------------------------------------------------------------------

% In order to perform the DIC analysis we use the NCC approach

tic

for k = 1:num_steps % Start at 1 until num_steps
    img1_ref = ImageSeries{k}; % The first image is the reference
    img2_def = ImageSeries{k + 1};  % The second image is the deformed

    j_node = 0; % Start the index counter

    for j = y_min:node_step_y:y_max % Loop over the Y-coordinates of the grid nodes
        j_node = j_node + 1; % Increment Y-index for the output array

        i_node = 0; % Reset X-index for each new row

        

        for i = x_min:node_step_x:x_max % Loop over the X-coordinates of the grid nodes
            i_node = i_node +1; % Increment X-index for the output array


            CC_max = -1; % Start max correlation for this node at minimum value

            % Extract the reference image motif, centered at (i, j)
            f_motif = img1_ref(j - corr_half:j + corr_half, i - corr_half : i + corr_half); % (row_start:row_end, col_start:col_end)

            % Calculate the constant reference image denominator (constant)

            f_denom_sq = sum(f_motif(:).^2); % computes the sum of squares of all values inside the array

            % Loop over the search range
            for u = -search_range_px:search_range_px
                for v = -search_range_px:search_range_px

                    % Extract the deformed image motif g, centered at (i+u,
                    % j+v)
                    g_motif = img2_def(j+v - corr_half:j+v + corr_half, i+u - corr_half : i+u + corr_half);


                    % ------- CALCULATE NCC --------

                    % Numerator: sum(f * g)
                    numerator = sum(sum(f_motif .* g_motif)); % matrix inner product of the two motifs

                    % Denominator: sqrt(Sum(f^2) * Sum(g^2))
                    g_denom_sq = sum(g_motif(:).^2);
                    denominator = sqrt(f_denom_sq * g_denom_sq); % total denominator

                    CC_now = numerator/denominator;

                    % We check for the maximum correlation
                    if CC_now > CC_max
                        CC_max = CC_now;
                        U_int_disp(j_node, i_node, k) = u; % Store the best integer u
                        V_int_disp(j_node, i_node, k) = v; % Store the best integer v
                        CC_max_store(j_node, i_node, k) = CC_max;
                    end
                end % end v search
            end    % end u search
        end        % end i-node
    end    % end j-node
end % end outermost loop

toc

%% 
%-------------------------------------------------------------------------
% Plotting displacement matrices
%-------------------------------------------------------------------------

% 1. Calculate the magnitude for step k
D_magnitude = sqrt(U_int_disp(:,:,k).^2 + V_int_disp(:,:,k).^2);

% 2. Plot the magnitude field
figure;
pcolor(X, Y, D_magnitude); % Use your X, Y grid coordinates
shading interp;             % Smooth the color transition
colormap jet;               % Use a good color map
colorbar;                   % Show the color scale
axis equal tight;           % Fit to the data
set(gca, 'YDir', 'normal'); % Match image orientation
title(['Displacement Magnitude (D) - Step ', num2str(k)]);

%% 
%-------------------------------------------------------------------------
% Calculating and plotting strain fields
%-------------------------------------------------------------------------

% Get the total number of steps to loop through
num_steps = size(U_int_disp, 3); 

% Define pre-allocation matrices

Epsilon_xx = zeros(NUM_nodes_y, NUM_nodes_x, num_steps); 
Epsilon_yy = zeros(NUM_nodes_y, NUM_nodes_x, num_steps); 
Gamma_xy = zeros(NUM_nodes_y, NUM_nodes_x, num_steps);

% Loop through each time step
for k = 1:num_steps
    % 1. Extract the displacement fields for the current step (k)
    U_k = U_int_disp(:,:,k); % U is displacement in X-direction
    V_k = V_int_disp(:,:,k); % V is displacement in Y-direction

    % 2. Calculate the Gradients using the node step sizes
    
    % [dU_dX, dU_dY] = gradient(U_matrix, spacing_X, spacing_Y)
    [dU_dX, dU_dY] = gradient(U_k, node_step_x, node_step_y);
    
    % [dV_dX, dV_dY] = gradient(V_matrix, spacing_X, spacing_Y)
    [dV_dX, dV_dY] = gradient(V_k, node_step_x, node_step_y);
    
    % ... proceed to strain calculation (Step 2)

    % Normal Strain Components
    Epsilon_xx(:,:,k) = dU_dX;
    Epsilon_yy(:,:,k) = dV_dY;

    % Engineering Shear Strain Component (Gamma_xy = dU/dY + dV/dX)
    Gamma_xy(:,:,k) = dU_dY + dV_dX;
end


% Define the strain component you want to animate (e.g., Epsilon_xx)
Strain_Component = Epsilon_xx; 
num_steps = size(Strain_Component, 3); 

figure; % Open a single figure window for the animation

% Set up the common plot properties once (axis limits, color map, etc.)
pcolor(X, Y, Strain_Component(:,:,1)); % Plot the first frame to initialize the plot object
shading interp;                       % Use smooth shading
colormap jet;                         % Use the jet color map (or 'parula' for a modern look)
colorbar;                             % Display the color scale
axis equal tight;                     % Ensure axes fit data nicely
set(gca, 'YDir', 'normal');           % Set image orientation

% Determine a consistent color scale (CLim) for all frames 
% This is CRITICAL so strain magnitude changes are visible, not just color shifts.
max_strain = max(Strain_Component(:));
min_strain = min(Strain_Component(:));
clim([min_strain max_strain]); % Lock the color limits across the animation

% Start the animation loop
fprintf('Starting strain animation...\n');
for k_plot = 1:num_steps
    % 1. Update the data in the existing plot object
    % Use 'CData' to quickly change the color data of the pcolor plot
    set(findobj(gca, 'Type', 'Surface'), 'CData', Strain_Component(:,:,k_plot));

    % 2. Update the title for the current frame
    title(['Normal Strain \epsilon_{xx}: Step ', num2str(k_plot), ...
           ' (Img ', num2str(k_plot), ' \rightarrow Img ', num2str(k_plot+1), ')']);
    
    % 3. Force MATLAB to draw the updated plot immediately
    drawnow;
    
    % Optional: Add a small pause to control the animation speed
     pause(0.5); 
end
fprintf('Animation complete.\n');
