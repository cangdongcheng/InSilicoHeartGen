% Robust Batch script to process biventricular meshes with error logging
clear;
clc;

%% 1. Configure paths
current_path = pwd; 
addpath(genpath(fullfile(current_path, '..', 'functions')));
addpath(genpath(fullfile(current_path, '..', 'dependencies')));

base_dir = 'D:\FYP_data\Yilin_20251222';

% Initialize error tracking
error_cases = {}; 

%% 2. Loop through all 25 samples (000 to 024)
for i = 0:24
    suffix = sprintf('%03d', i);
    folder_name = ['Average_merge_', suffix];
    
    working_dir = fullfile(base_dir, folder_name);
    name_origin = [folder_name, '.vtk'];
    name_final = ['OUT_', folder_name];
    case_folder = fullfile(working_dir, name_final);
    
    fprintf('\n==================================================\n');
    fprintf('Checking Sample %d/24: %s\n', i, folder_name);
    
    if exist(case_folder, 'dir')
        fprintf('Output folder already exists. Skipping %s...\n', folder_name);
        continue;
    end
    
    if ~exist(working_dir, 'dir')
        warning('Directory not found: %s. Skipping...', working_dir);
        continue;
    end
    
    % Use try-catch to handle TetGen or MATLAB crashes gracefully
    try
        % Create folder only when we start processing
        mkdir(case_folder);
        
        %% 3. Mesh Generation Pipeline
        input_file = fullfile(working_dir, name_origin);
        if ~exist(input_file, 'file')
            error('Surface mesh not found.');
        end
        
        disp(['Reading surface mesh: ', name_origin]);
        surf0 = vtkRead(input_file);
        
        % --- MESH PRE-REPAIR: Cleaning duplicate points ---
        % This can often prevent the 'Assertion failed' TetGen error
        [unique_pts, ~, ic] = unique(surf0.points, 'rows');
        surf0.points = unique_pts;
        surf0.cells = ic(surf0.cells); 
                   
        %% Unit Scaling Check
        dist = sqrt(sum((surf0.points - repmat(surf0.points(1,:), length(surf0.points), 1)).^2, 2));
        scaledist = max(dist);
        
        if scaledist > 5 && scaledist <= 50
            surf0.points = surf0.points .* 10;
        end

        %% Generate Coarse (1.5)
        disp('Generating Coarse Tetrahedral Mesh...');
        MeshCoarse = tetrahedral_meshing(surf0, 1.5, [], []);
        MeshCoarse.points = MeshCoarse.points ./ 10;
        vtkWrite(MeshCoarse, fullfile(case_folder, 'Coarse.vtu'));

        %% Generate Fine (1.0)
        disp('Generating Fine Tetrahedral Mesh...');
        MeshFine = tetrahedral_meshing(surf0, 1.0, [], []);
        MeshFine.points = MeshFine.points ./ 10;
        vtkWrite(MeshFine, fullfile(case_folder, 'fine.vtu'));

        %% Flip Normals
        disp('Verifying volumetric orientation...');
        sur_coarse = vtkDataSetSurfaceFilter(MeshCoarse);
        TR_Surf = triangulation(double(sur_coarse.cells), double(sur_coarse.points));
        TR_Vol = triangulation(double(MeshCoarse.cells), double(MeshCoarse.points));
        centroid = meshcentroid(double(sur_coarse.points), double(sur_coarse.cells));
        facenormals = faceNormal(TR_Surf);
        offsetPoint = centroid(1,:) + 1e-4 * facenormals(1,:);
        
        if ~isnan(pointLocation(TR_Vol, offsetPoint))
            MeshCoarse.cells(:,[3 4]) = MeshCoarse.cells(:,[4 3]);
            MeshFine.cells(:,[3 4]) = MeshFine.cells(:,[4 3]);
            vtkWrite(MeshCoarse, fullfile(case_folder, 'Coarse.vtu'));
            vtkWrite(MeshFine, fullfile(case_folder, 'fine.vtu'));
        end
        
        disp(['--- SUCCESS: ', folder_name, ' ---']);

    catch ME
        % If an error occurs, log it and clean up
        fprintf('!! ERROR encountered in %s: %s\n', folder_name, ME.message);
        error_cases{end+1} = folder_name; %#ok<AGROW>
        
        % Cleanup: Remove the incomplete directory so it can be re-tried later
        if exist(case_folder, 'dir')
            rmdir(case_folder, 's');
        end
        continue; % Move to next iteration
    end
end

%% 4. Final Summary
fprintf('\n==================================================\n');
fprintf('BATCH PROCESSING COMPLETE\n');
if isempty(error_cases)
    fprintf('All cases processed successfully.\n');
else
    fprintf('The following %d cases FAILED:\n', length(error_cases));
    for k = 1:length(error_cases)
        fprintf('- %s\n', error_cases{k});
    end
end
fprintf('==================================================\n');