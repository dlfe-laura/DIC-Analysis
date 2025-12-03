%% 
%-------------------------------------------------------------------------
% First we define the parameters
%-------------------------------------------------------------------------

    % Node points. We define a rectangular grid with limits nd step sizes
    x_min = 61;
    x_max = 1056;
    node_step_x = 10; % the step size between analysis node

    y_min = 601;
    y_max = 1203;
    node_step_y = 10;
    
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
    search_range_px = 3;
    
    % Correlation window size (the specific size of a pattern I am looking for). 
    % The template will slide across my entire analysis window to find matches.
    corr_half = 12;
    ref_x= 21 ;
    ref_y = 21 ;
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

% Define options for fminsearch once (increases limits to prevent early exit)
options = optimset('Display','off', 'MaxFunEvals', 1000, 'MaxIter', 1000);


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

            % We need to store the correlation surface for this node (pre-allocation subpixel)
            search_size = 2 * search_range_px + 1;
            CCmat = zeros(search_size, search_size);
            u_index = 0;

            % Loop over the search range
            for u = -search_range_px:search_range_px
                u_index = u_index + 1; 
                v_index = 0; % Reset v_index for every new column (u)

                for v = -search_range_px:search_range_px
                    v_index = v_index + 1; % Increment v_index for every row (v)

                    % Extract the deformed image motif g, centered at (i+u,j+v)
                    g_motif = img2_def(j+v - corr_half:j+v + corr_half, i+u - corr_half : i+u + corr_half);


                    % ------- CALCULATE NCC --------

                    % Numerator: sum(f * g)
                    numerator = sum(sum(f_motif .* g_motif)); % matrix inner product of the two motifs

                    % Denominator: sqrt(Sum(f^2) * Sum(g^2))
                    g_denom_sq = sum(g_motif(:).^2);
                    denominator = sqrt(f_denom_sq * g_denom_sq); % total denominator

                    CC_now = numerator/denominator;

                    % Store in Matrix for Sub-pixel fit
                    CCmat(v_index, u_index) = CC_now;

                    % We check for the maximum correlation
                    if CC_now > CC_max
                        CC_max = CC_now;
                        U_int_disp(j_node, i_node, k) = u; % Store the best integer u
                        V_int_disp(j_node, i_node, k) = v; % Store the best integer v
                        CC_max_store(j_node, i_node, k) = CC_max;
                    end
                end % end v search
            end    % end u search

            % ============================================================
            %  SUB-PIXEL REFINEMENT
            % ============================================================
            try
                % 1. Create the grid for interpolation
                % Transpose inputs to match MATLAB matrix standard
                [SXI, SYI] = meshgrid(-search_range_px:search_range_px, -search_range_px:search_range_px);
                
                % 2. Create the Interpolant
                % Using 'pchip' or 'makima' is often safer than 'spline' to avoid overshooting, 
                % but we will stick to 'spline' with a bounds check.
                CCgrid = griddedInterpolant(SXI', SYI', CCmat', 'spline');
                
                % 3. Define the start point (Best Integer Result)
                start_point = [U_int_disp(j_node, i_node, k), V_int_disp(j_node, i_node, k)];
                
                % 4. Run fminsearch
                [xy_sub, fval] = fminsearch(@(xy) -CCgrid(xy(1), xy(2)), start_point, options);
                
                % --- SANITY CHECK (The Fix) ---
                % If the optimizer found a point outside our search window (-10 to 10),
                % it is a mathematical artifact. Discard it.
                u_refined = xy_sub(1);
                v_refined = xy_sub(2);
                
                if abs(u_refined) <= search_range_px && abs(v_refined) <= search_range_px
                    % Result is valid, store it
                    U_int_disp(j_node, i_node, k) = u_refined;
                    V_int_disp(j_node, i_node, k) = v_refined;
                else
                    % Result is Garbage (10^69), keep the Integer result!
                    % Do nothing here.
                end
                
            catch
                % If fminsearch crashes, keep Integer result.
            end
            % ============================================================
        end        % end i-node
    end    % end j-node
end % end outermost loop

toc

%% 
%-------------------------------------------------------------------------
% Plotting Displacement matrices
%-------------------------------------------------------------------------
fprintf('Preparing displacement animation...\n');

% 1. Pre-calculate Magnitude for ALL steps to establish a global color scale
% This creates a 3D matrix where the 3rd dimension is time
D_total = sqrt(U_int_disp.^2 + V_int_disp.^2);

% Find the maximum displacement across the entire series for the color bar
max_disp = max(D_total(:));
min_disp = min(D_total(:));

% 2. Setup the Figure (Initialize with the first frame)
figure;
hDispPlot = pcolor(X, Y, D_total(:,:,1)); % Save the plot handle 'hDispPlot'
shading interp;             % Smooth colors
colormap jet;               % Rainbow map
colorbar;                   % Show scale
axis equal tight;           % Fit to data
set(gca, 'YDir', 'normal'); % Correct orientation

% Lock the color limits so colors represent consistent values over time
if max_disp > min_disp
    clim([min_disp, max_disp]); 
end

xlabel('X-coordinate (pixels)');
ylabel('Y-coordinate (pixels)');

% 3. Run the Animation Loop
num_steps = size(U_int_disp, 3);

fprintf('Starting displacement animation...\n');
for k = 1:num_steps
    % --- Safety Check: Stop if window is closed ---
    if ~isvalid(hDispPlot)
        break; 
    end
    
    % Update the plot data using the pre-calculated magnitude
    set(hDispPlot, 'CData', D_total(:,:,k));
    
    % Update the title
    title(['Displacement Magnitude (D) - Step ', num2str(k), ...
           ' (Sub-pixel Refined)']);
    
    % Render the frame
    drawnow;
    
    % Pause briefly to control playback speed (adjust as needed)
    pause(0.5); 
end
fprintf('Displacement animation complete.\n');

%% 
%-------------------------------------------------------------------------
% Calculating strain fields 
%-------------------------------------------------------------------------
fprintf('Calculating strain fields...\n');

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
    
    % Calculate the Gradients using the node step sizes
    % [dU/dx, dU/dy]
    [dU_dX, dU_dY] = gradient(U_k, node_step_x, node_step_y);
    
    % [dV/dx, dV/dy]
    [dV_dX, dV_dY] = gradient(V_k, node_step_x, node_step_y);
    
    % Calculate Small Strain Tensor Components
    % Normal Strain X (E_xx = dU/dx)
    Epsilon_xx(:,:,k) = dU_dX;
    
    % Normal Strain Y (E_yy = dV/dy)
    Epsilon_yy(:,:,k) = dV_dY;
    
    % Engineering Shear Strain (Gamma_xy = dU/dy + dV/dx)
    Gamma_xy(:,:,k) = dU_dY + dV_dX;
end

fprintf('Strain calculation complete.\n');

%% 
%-------------------------------------------------------------------------
% Plotting strain fields (Animation)
%-------------------------------------------------------------------------

% Change this variable to plot different strains: Epsilon_xx, Epsilon_yy, or Gamma_xy
Strain_Component = Epsilon_yy; % Select component to plot
Strain_Name = 'Normal Strain \epsilon_{xx}'; % Title for the plot (change)

num_steps = size(Strain_Component, 3); 

figure; % Open a new figure window
% Setup the plot with the first frame
hPlot = pcolor(X, Y, Strain_Component(:,:,1)); 
shading interp;             % Interpolates colors for a smooth look
colormap jet;               % Jet is standard, try 'parula' or 'turbo' for better contrast
colorbar;                   
axis equal tight;           % Fit axes to data
set(gca, 'YDir', 'normal'); % Correct image orientation

xlabel('X (pixels)'); % Labelling axis
ylabel('Y (pixels)');

% 2. ROBUST COLOR SCALING
% Sub-pixel noise can sometimes create a single "hot" pixel that ruins the scale.
% We calculate the global max/min, but we also check for outliers.

% Find global max and min across all time steps
global_max = max(Strain_Component(:));
global_min = min(Strain_Component(:));

% Calculate a symmetric limit
c_limit = max(abs(global_min), abs(global_max));

clim([-c_limit, c_limit]); % Set symmetric limits centered on zero

% Animation loop
fprintf('Starting strain animation...\n');

for k_plot = 1:num_steps
    % Update the color data (CData) of the existing plot object
    set(hPlot, 'CData', Strain_Component(:,:,k_plot));
    
    % Update title
    title([Strain_Name, ': Step ', num2str(k_plot), ...
           ' (Img ', num2str(k_plot), ' \rightarrow Img ', num2str(k_plot+1), ')']);
    
    % Update the drawing
    drawnow;
    
    % Adjust speed (seconds)
    pause(0.4); 
end
fprintf('Animation complete.\n');
