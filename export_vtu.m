% export_vtu.m
% Builds a VTU from CSV files and writes it to the same folder.
% Usage: set ensi_path below, then run.
%   ensi_path = 'C:\...\outputs\NUH_test\1\ensi1'

clear
current_path = pwd;
addpath(genpath(fullfile(current_path, 'functions')));
addpath(genpath(fullfile(current_path, 'dependencies')));

%% ---- SET THIS PATH ----
ensi_path = 'C:\Users\e1032484\InSilicoHeartGen\outputs\NUH_test\1\ensi_Fine_1';
%% ------------------------

csv_dir = fullfile(ensi_path, 'CSVFiles');
assert(isfolder(csv_dir), 'CSVFiles folder not found in: %s', ensi_path);

% Find the prefix (e.g. "1_coarse") from the xyz file
xyz_file = dir(fullfile(csv_dir, '*_xyz.csv'));
assert(~isempty(xyz_file), 'No *_xyz.csv found in %s', csv_dir);
prefix = erase(xyz_file(1).name, '_xyz.csv');
fprintf('Prefix: %s\n', prefix);

%% Build mesh from CSVs
points = csvread(fullfile(csv_dir, [prefix '_xyz.csv']));
cells  = csvread(fullfile(csv_dir, [prefix '_tetra.csv']));
n_nodes = size(points, 1);
n_cells = size(cells, 1);
fprintf('Mesh: %d nodes, %d cells\n', n_nodes, n_cells);

mesh.points = points;
mesh.cells  = int32(cells);
mesh.cellTypes = repmat(uint8(10), n_cells, 1);  % VTK_TETRA

%% Point data (per-node fields)
node_fields = {
    'Fiber',           'nodefield_fibre'
    'Sheet',           'nodefield_sheet'
    'Normal',          'nodefield_normal'
    'Cobiveco',        'cobiveco'
    'Apex2Base',       'nodefield_cobiveco-ab'
    'Transmurality',   'nodefield_cobiveco-tm'
    'Rotational',      'nodefield_cobiveco-rt'
    'Transventricular','nodefield_cobiveco-tv'
    'RV2LV_geodesic',  'nodefield_cobiveco-rvlv'
    'RV2LV_projection','nodefield_cobiveco-rvlv-projection'
    'Ant2Post',        'nodefield_cobiveco-aprt'
};

for i = 1:size(node_fields, 1)
    fname = fullfile(csv_dir, [prefix '_' node_fields{i,2} '.csv']);
    if isfile(fname)
        val = csvread(fname);
        if size(val, 1) == n_nodes
            mesh.pointData.(node_fields{i,1}) = double(val);
            fprintf('  + pointData.%s [%dx%d]\n', node_fields{i,1}, size(val));
        else
            warning('Skipping pointData.%s: %d rows ~= %d nodes', node_fields{i,1}, size(val,1), n_nodes);
        end
    else
        fprintf('  - pointData.%s not found, skipping\n', node_fields{i,1});
    end
end

% Boundary node fields → binary per-node arrays
boundary_fields = {
    'EP_LVnodes', 'boundarynodefield_ep-lvnodes'
    'EP_RVnodes', 'boundarynodefield_ep-rvnodes'
};

for i = 1:size(boundary_fields, 1)
    fname = fullfile(csv_dir, [prefix '_' boundary_fields{i,2} '.csv']);
    if isfile(fname)
        idx = csvread(fname);
        flag = zeros(n_nodes, 1);
        valid = idx(idx >= 1 & idx <= n_nodes);
        flag(valid) = 1;
        mesh.pointData.(boundary_fields{i,1}) = double(flag);
        fprintf('  + pointData.%s [%d flagged nodes]\n', boundary_fields{i,1}, numel(valid));
    end
end

%% Cell data (per-element fields)
cell_fields = {
    'Material', 'material_tetra'
};

for i = 1:size(cell_fields, 1)
    fname = fullfile(csv_dir, [prefix '_' cell_fields{i,2} '.csv']);
    if isfile(fname)
        val = csvread(fname);
        if size(val, 1) == n_cells
            mesh.cellData.(cell_fields{i,1}) = double(val);
            fprintf('  + cellData.%s [%dx%d]\n', cell_fields{i,1}, size(val));
        else
            warning('Skipping cellData.%s: %d rows ~= %d cells', cell_fields{i,1}, size(val,1), n_cells);
        end
    end
end

%% Write VTU
out_file = fullfile(ensi_path, [prefix '_fields.vtu']);
vtkWrite(mesh, out_file);
fprintf('Written: %s\n', out_file);
