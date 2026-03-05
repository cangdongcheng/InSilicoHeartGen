clear;
clc;

%% 1. Configure paths
current_path = pwd; 
addpath(genpath(fullfile(current_path, 'functions')));
addpath(genpath(fullfile(current_path, 'dependencies')));

root_dir = 'D:\VTK_Merged_ED';
log_file = fullfile(root_dir, 'Processing_Log.txt');
dir_content = dir(root_dir);

% Filter for valid directories
valid_folders = dir_content([dir_content.isdir] & ~ismember({dir_content.name}, {'.', '..'}));
total_cases = length(valid_folders);
failed_cases = {}; 

%% 2. Loop through each item
for i = 1:total_cases
    case_name = valid_folders(i).name;
    working_dir = fullfile(root_dir, case_name);
    
    % Update progress line in Command Window
    percent = (i / total_cases) * 100;
    fprintf('\rProgress: [%3.0f%%] - Processing: %s', percent, case_name);
    
    %% 3. Define filenames
    name_origin = [case_name, '.vtk']; 
    case_folder = fullfile(working_dir, 'Volume_OUT');
    
    if ~exist(case_folder, 'dir')
        mkdir(case_folder);       
    end
    
    %% 4. Mesh Generation Pipeline
    input_file = fullfile(working_dir, name_origin);
    
    if ~exist(input_file, 'file')
        failed_cases{end+1} = sprintf('%s: .vtk file not found', case_name);
        continue;
    end
    
    try
        % Execute meshing silently
        evalc([...
            'surf0 = vtkRead(input_file);' ...
            'dist = sqrt(sum((surf0.points - repmat(surf0.points(1,:), size(surf0.points,1), 1)).^2, 2));' ...
            'if max(dist) > 5 && max(dist) <= 50, surf0.points = surf0.points .* 10; end;' ...
            'MeshCoarse = tetrahedral_meshing(surf0, 1.5, [], []);' ...
            'MeshCoarse.points = MeshCoarse.points ./ 10;' ...
            'coarse_path = fullfile(case_folder, ''Coarse.vtu'');' ...
            'vtkWrite(MeshCoarse, coarse_path);' ...
            ...
            'sur_coarse = vtkDataSetSurfaceFilter(MeshCoarse);' ...
            'TR_Surf = triangulation(double(sur_coarse.cells), double(sur_coarse.points));' ...
            'TR_Vol = triangulation(double(MeshCoarse.cells), double(MeshCoarse.points));' ...
            'centroid = meshcentroid(double(sur_coarse.points), double(sur_coarse.cells));' ...
            'facenormals = faceNormal(TR_Surf);' ...
            'offsetPoint = centroid(1,:) + 1e-4 * facenormals(1,:);' ...
            'if ~isnan(pointLocation(TR_Vol, offsetPoint)),' ...
            'MeshCoarse.cells(:,[3 4]) = MeshCoarse.cells(:,[4 3]);' ...
            'vtkWrite(MeshCoarse, coarse_path); end;' ...
        ]);
        
    catch ME
        % Clean the error message for the log
        clean_msg = regexprep(ME.message, '[\n\r]+', ' ');
        failed_cases{end+1} = sprintf('%s: %s', case_name, clean_msg);
    end
end

%% 5. Write Log File and Final Report
fid = fopen(log_file, 'w');
fprintf(fid, 'Processing Log - %s\n', datestr(now));
fprintf(fid, 'Total cases found: %d\n', total_cases);
fprintf(fid, '-------------------------------------------\n');

if isempty(failed_cases)
    fprintf(fid, 'Status: All cases processed successfully.\n');
    fprintf('\n\nProcessing Complete. All cases successful.\n');
else
    fprintf(fid, 'Status: %d cases failed.\n\nFAILED CASES:\n', length(failed_cases));
    for k = 1:length(failed_cases)
        fprintf(fid, '- %s\n', failed_cases{k});
    end
    fprintf('\n\nProcessing Complete. %d errors logged to %s\n', length(failed_cases), log_file);
end
fclose(fid);