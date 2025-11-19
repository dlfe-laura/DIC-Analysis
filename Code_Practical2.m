%-------------------------------------------------------------------------
% First we define the parameters
    % Node points. We define a rectangular grid with limits nd step sizes

    x_min = 0;
    x_max = 822;
    x_step = 10; % the step size between analysis node

    y_min = 0;
    y_max = 548;
    y_step = 10;
    
    % We define the vectors for each direction to define the 2D-matrices

    x_vect = x_min:x_step:x_max; % = start:step:end
    y_vect = y_min:y_step:y_max;
    
    %The 2D matrices

    [X, Y] = meshgrid(x_vect, y_vect);

    % 4. (Optional) Visualize the grid nodes
    %figure;
    %plot(X(:), Y(:), 'b.'); % Plot all grid points as blue dots
    %axis equal;             % Ensures dx and dy steps look proportional
    %xlabel('X-coordinate');
    %ylabel('Y-coordinate');
    %title('Regular Rectangular Grid');
    %grid on;

    % Search window size (smaller window that moves around inside my
    % analysis window to perform a local calculation at each grid node.)

    x_wd_min = 0;
    x_wd_max = 411; % x matrix divided by 2

    y_xd_min = 0;
    y_xd_max = 274; % y matrix divided by 2
    
    % Correlation window size (the specific size of a pattern I am
    % looking for). The template will slide across my entire 
    % Analysis Window to find matches.

    ref_x= 32 ; % limit of the grid. maybe a bit bigger (32x32) but take 
    % into account later the resizing
    ref_y = 32 ;

    % I have to define the correlation coefficient (?) CHECK THIS

%-------------------------------------------------------------------------
% Read, plot image, check how many pixels it has, adjust my parameters and
% limits and then resize the images. After start looping

    % Read test image to get information for parameters

    myFolder = 'C:\Users\la3314de\Downloads\drive-download-20251113T090234Z-1-001'; 

    fileName = 'success rl-0000_0.tif';

    try
    % 2. Read the Image
    % imread() loads the image data into a matrix
    img = imread(fileName);
    
    % 3. Plot the Image
    %figure;      % Create a new figure window
    %imshow(img); % Display the image
    %title(['Image: ', fileName]); % Add a title
    
    % 4. Measure the Size
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
            fprintf('Created output folder: %s\n', outputFolder);
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

        fprintf('Found %d TIF images. Starting resize...\n', numFiles);

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
                % imresize(image, scale)
                resizedImg = imresize(img, scaleFactor);
        
                % Get new [height, width]
                newSize = [size(resizedImg, 1), size(resizedImg, 2)];
        
                % --- Save the new image ---
                imwrite(resizedImg, newFullFileName);
        
                % --- Log the results ---
                % "Effective Pixel Size Change"
                % If scaleFactor = 0.25 (smaller), one new pixel
                % covers the area of 2x2 old pixels.
                % The effective size of the new pixel is 1 / 0.25 = 4.
                pixelChange = 1 / scaleFactor;
        
                resultsTable(k, :) = {baseFileName, ...
                                    mat2str(origSize), ...
                                    mat2str(newSize), ...
                                    pixelChange};
                              
                fprintf('Processed: %s\n', baseFileName);
        
            catch ME
                fprintf('Error processing %s: %s\n', baseFileName, ME.message);
                resultsTable(k, :) = {baseFileName, 'ERROR', 'ERROR', NaN};
            end
        end

        % Display the results of the resizing

        fprintf('\n--- Batch Resize Complete ---\n');
        disp('Results:');
        disp(resultsTable);


%-------------------------------------------------------------------------
% Loop over DIC analysis nodes (grid) 