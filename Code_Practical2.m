%-------------------------------------------------------------------------
% First we define the parameters
%-------------------------------------------------------------------------

    % Node points. We define a rectangular grid with limits nd step sizes
    x_min = 41;
    x_max = 1056;
    node_step_x = 10; % the step size between analysis node

    y_min = 41;
    y_max = 1604;
    node_step_y = 10;
    
    % We define the vectors for each direction to define the 2D-matrices
    x_vect = x_min:node_step_x:x_max; % = start:step:end
    y_vect = y_min:node_step_y:y_max;
    
    % The number of nodes
    NUM_nodes_x = length(x_vect);
    NUM_nodes_y = length(y_vect);

    %The 2D matrices
    [X, Y] = meshgrid(x_vect, y_vect);

    % Search window size (smaller window that moves around inside my
    % analysis window to perform a local calculation at each grid node.)
    % It defines the range of possible integer displacements (u, v) around
    % node x_i
    search_range_px = 20;
    
    % Correlation window size (the specific size of a pattern I am
    % looking for). The template will slide across my entire 
    % Analysis Window to find matches.
    corr_half = 20;
    ref_x= 41 ;
    ref_y = 41 ;

    % Prealocate matrices to store the results
    U_int_disp = zeros(NUM_nodes_y, NUM_nodes_x); % Integer displacement in X
    V_int_disp = zeros(NUM_nodes_y, NUM_nodes_x); % Integer displacement in Y
    CC_max_store = zeros(NUM_nodes_y, NUM_nodes_x); % Max correlation coefficient

%-------------------------------------------------------------------------
% Read, plot image, check how many pixels it has, adjust my parameters and
% limits and then resize the images. After start looping
%-------------------------------------------------------------------------

    % Read test image to get information for parameters

    myFolder = 'C:\Users\la3314de\Downloads\drive-download-20251113T090234Z-1-001'; 

    fileName = 'success rl-0000_0.tif';

    try
    % Read the Image
    % imread() loads the image data into a matrix
    img = imread(fileName);
    
    % Measure the Size
    % size() returns the dimensions of the matrix
    [height, width] = size(img);
    
    % Display the size in the Command Window
    disp(['--- Image Size (', fileName, ') ---']);
    disp(['Height: ', num2str(height), ' pixels']);
    disp(['Width:  ', num2str(width), ' pixels']);

    catch ME
        % Handle errors, like "file not found"
        fprintf('Error reading image "%s":\n', fileName);
        fprintf('%s\n', ME.message);
        fprintf('Please check that the file exists and is in the correct path.\n');
    end

    % Resize my images for better handling
    scaleFactor = 0.25;
        % Create a new folder for the resized images
        outputFolder = fullfile(myFolder, 'resized_images');
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
            fprintf('Created output folder: %s\n', outputFolder); % This is okay to keep for setup info
        end
        % Get a list of all .tif files
        filePattern = fullfile(myFolder, '*.tif');
        imageList = dir(filePattern);
        % Check if files were found
        if isempty(imageList)
            warning('No TIF files found in this folder: %s', myFolder);
            return;
        end
        % Pre-allocate a table to store the results
        numFiles = length(imageList);
        resultsTable = table('Size', [numFiles, 4], ...
            'VariableTypes', {'string', 'string', 'string', 'double'}, ...
            'VariableNames', {'FileName', 'OriginalSize', 'NewSize', 'PixelSizeChangeFactor'});
        
        % REMOVED: fprintf('Found %d TIF images. Starting resize...\n', numFiles);

        % Loop, resize and save the new images
        for k = 1:numFiles
            % --- Get file paths ---
            baseFileName = imageList(k).name;
            fullFileName = fullfile(myFolder, baseFileName);
            newFullFileName = fullfile(outputFolder, baseFileName);
    
            try
                % --- Read and get original size ---
                img = imread(fullFileName);
                % Get [height, width]
                origSize = [size(img, 1), size(img, 2)];
        
                % --- Resize the image ---
                resizedImg = imresize(img, scaleFactor);
        
                % Get new [height, width]
                newSize = [size(resizedImg, 1), size(resizedImg, 2)];
        
                % --- Save the new image ---
                imwrite(resizedImg, newFullFileName);
        
                % --- Log the results (Silently store) ---
                pixelChange = 1 / scaleFactor;
        
                resultsTable(k, :) = {baseFileName, ...
                                    mat2str(origSize), ...
                                    mat2str(newSize), ...
                                    pixelChange};
                              
                % REMOVED: fprintf('Processed: %s\n', baseFileName);
        
            catch ME
                % KEPT: Only print when an error occurs
                fprintf('Error processing %s: %s\n', baseFileName, ME.message);
                resultsTable(k, :) = {baseFileName, 'ERROR', 'ERROR', NaN};
            end
        end
        
        % REMOVED: Block that displays the final results table
        % fprintf('\n--- Batch Resize Complete ---\n');
        % disp('Results:');
        % disp(resultsTable);


        % --- Read the necessary images for DIC ---
        % You need to know the specific filenames of the reference (img1) and deformed (img2) images.
        
        ref_filename = 'success rl-0000_0.tif'; % Example reference file
        def_filename = 'success rl-0001_0.tif'; % Example deformed file

        % Load the images from the 'resized_images' folder
        img1_ref = imread(fullfile(outputFolder, ref_filename)); 
        img2_def = imread(fullfile(outputFolder, def_filename));

        % Convert images to double precision (REQUIRED for accurate NCC/ZNCC calculations)
        img1_ref = double(img1_ref);
        img2_def = double(img2_def);

%-------------------------------------------------------------------------
% Loop over DIC analysis nodes (grid) 
%-------------------------------------------------------------------------

% In order to perform the DIC analysis we use the NCC approach

j_node = 0; % Start the index counter

for j = y_min:node_step_y:y_max % Loop over the Y-coordinates of the grid nodes
    j_node = j_node + 1; % Increment Y-index for the output array

    x_node = 0; % Reset X-index for each new row

    for i = x_min:node_step_x:x_max % Loop over the X-coordinates of the grid nodes
        x_node = x_node +1; % Increment X-index for the output array


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
                    U_int_disp(j_node, i_node) = u; % Store the best integer u
                    V_int_disp(j_node, i_node) = v; % Store the best integer v
                    CC_max_store(j_node, i_node) = CC_max;
                end
            end % end v search
        end    % end u search
    end        % end i (x-node)      
end    % end j (y-node)

%-------------------------------------------------------------------------
% Plotting
%-------------------------------------------------------------------------

figure;
surf(X, Y, U_int_disp, 'EdgeColor', 'none'); 
view(2); % Set the view to 2D (top-down) for a clear map
colorbar;
axis equal tight;
xlabel('X-coordinate (pixels)');
ylabel('Y-coordinate (pixels)');
title('Surface Map of Horizontal Displacement (U)');
