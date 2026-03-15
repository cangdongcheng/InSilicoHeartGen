% This is script which runs the pipeline from an existing biventricular surface mesh to simulation files
clear;
clc;

%% 1. Configure paths
current_path = pwd; % Should be C:\Users\e1032484\InSilicoHeartGen
addpath(genpath(fullfile(current_path, 'functions')));
addpath(genpath(fullfile(current_path, 'dependencies')));

% Define the working directory on D: drive
working_dir = 'D:\FYP_data\Yilin_20251222\Average_merge_008';

%% 2. Define filenames
name_origin = 'Average_merge_008.vtk'; % Updated extension to .vtk
name_final = 'OUT_Average_merge_008';

% Create the output subfolder inside the D: drive directory
case_folder = fullfile(working_dir, name_final);
if ~exist(case_folder, 'dir')
    mkdir(case_folder);       
end

%% 3. Mesh Generation Pipeline
for index = 1
    %% Read original surface mesh
    input_file = fullfile(working_dir, name_origin);
    if ~exist(input_file, 'file')
        error('Surface mesh not found at: %s', input_file);
    end
    
    disp(['Reading surface mesh: ', name_origin]);
    surf0 = vtkRead(input_file);        
               
    %% Unit Scaling Check
    % Check if mesh is in mm or cm and convert to mm for meshing
    dist = sqrt(sum((surf0.points - repmat(surf0.points(1,:), length(surf0.points), 1)).^2, 2));
    scaledist = max(dist); 
    
    if scaledist > 50 && scaledist < 500
        disp('Scale detection: mm');
    elseif scaledist > 5 && scaledist <= 50
        disp('Scale detection: cm -> Converting to mm for processing');
        surf0.points = surf0.points .* 10;
    else
        warning('Mesh dimensions are unusual. Proceeding with raw values.');
    end

    %% Generate Coarse Tetrahedral Volume Mesh
    disp('Generating Coarse Tetrahedral Mesh (1.5 resolution)...');
    mesh_resolution_coarse = 1.5;
    MeshCoarse = tetrahedral_meshing(surf0, mesh_resolution_coarse, [], []);
    
    % Scale back to cm for simulation compatibility
    MeshCoarse.points = MeshCoarse.points ./ 10;
    vtkWrite(MeshCoarse, fullfile(case_folder, 'Coarse.vtu'));

    %% Generate Fine Tetrahedral Volume Mesh
    disp('Generating Fine Tetrahedral Mesh (1.0 resolution)...');
    mesh_resolution_fine = 1.0;
    MeshFine = tetrahedral_meshing(surf0, mesh_resolution_fine, [], []);    
    MeshFine.points = MeshFine.points ./ 10;
    
    % Save fine mesh
    vtkWrite(MeshFine, fullfile(case_folder, 'fine.vtu'));

    %% Flip Normals (if necessary)
    disp('Verifying volumetric orientation...');
    sur_coarse = vtkDataSetSurfaceFilter(MeshCoarse);
    TR_Surf = triangulation(double(sur_coarse.cells), double(sur_coarse.points));
    TR_Vol = triangulation(double(MeshCoarse.cells), double(MeshCoarse.points));
    
    centroid = meshcentroid(double(sur_coarse.points), double(sur_coarse.cells));
    facenormals = faceNormal(TR_Surf);
    offsetPoint = centroid(1,:) + 1e-4 * facenormals(1,:);
    
    if ~isnan(pointLocation(TR_Vol, offsetPoint))
        disp('Flipping inward-facing normals...');
        MeshCoarse.cells(:,[3 4]) = MeshCoarse.cells(:,[4 3]);
        MeshFine.cells(:,[3 4]) = MeshFine.cells(:,[4 3]);
        vtkWrite(MeshCoarse, fullfile(case_folder, 'Coarse.vtu'));
        vtkWrite(MeshFine, fullfile(case_folder, 'fine.vtu'));
    end

    disp('--- SUCCESS: Volume meshes generated in OUT_Average_merge_008 folder ---');
end