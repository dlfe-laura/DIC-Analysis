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

x_min = 200;
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

num_frames = numFiles; % Count number of frames
num_steps = num_frames - 5; % Calculate the number of analysis steps (at step 30 rupture appears)

% We create zero matrices to store the displacement for each node
DX = zeros(NUM_nodes_y, NUM_nodes_x);
DY = zeros(NUM_nodes_y, NUM_nodes_x);

% Pre-allocate cell arrays for history
DX_all = cell(numFiles, 1);
DY_all = cell(numFiles, 1);

% Set Frame 1 to zero (since it is the reference)
DX_all{1} = DX; 
DY_all{1} = DY;

%%

%-------------------------------------------------------------------
% Loop for image correlation
%-------------------------------------------------------------------

tic

for k = 1:num_steps % Start at 1 until num_steps. It just loops through all the images

    img1_ref = ImageSeries{k}; % We define our reference image
    img2_def = ImageSeries{k+1}; % We define our deformed image

    dx_inc = zeros(NUM_nodes_y, NUM_nodes_x);     % Displacement in x
    dy_inc = zeros(NUM_nodes_y, NUM_nodes_x);     % Displacement in y
    CC_store = zeros(NUM_nodes_y, NUM_nodes_x);   % Zero matrix to store the CC

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

                    % Verifies that if x_min_d < 1 OR x_max_d > W ... it
                    % will skip that iteration of the loop
                    if x_min_d < 1 || x_max_d > W || y_min_d < 1 || y_max_d > H
                        continue;
                    end

                    zone_d = img2_def(y_min_d:y_max_d, x_min_d:x_max_d); % Zone of the deformed image

                    % Now we check the CC with the CC_max using NCC

                    % Numerator
                    numerator = sum(sum(zone_0 .* zone_d)); % Matrix inner product of the two images

                    % Denominator
                    denominator = sqrt (sum (sum (zone_0 .^2)) * sum (sum(zone_d .^ 2)));

                    % CC
                    CC = numerator/denominator;

                    % Check for maximum correlation
                    if CC > CC_max % If CC is greater than CC_max then store the value and its displacements
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
    
    % Store the new total in the cell array for Frame k
    DX_all{k} = DX;
    DY_all{k} = DY;

end % Image sequence loop
    
toc
