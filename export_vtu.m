% export_vtu.m
% Attaches all fields from Case_coarse.mat onto Coarse.vtu and writes a
% single VTU file openable in ParaView.

clear
current_path = pwd;
addpath(genpath(fullfile(current_path, 'functions')));
addpath(genpath(fullfile(current_path, 'dependencies')));

%% Configure
case_number = 1;
case_folder  = fullfile(current_path, 'outputs', num2str(case_number));
ensi_folder  = fullfile(case_folder, ['ensi', num2str(case_number)]);
mesh_file    = fullfile(case_folder, 'Coarse.vtu');
mat_file     = fullfile(ensi_folder, 'Case_coarse.mat');
out_file     = fullfile(case_folder, ['Case_coarse_fields_', num2str(case_number), '.vtu']);

%% Load verified mesh
mesh = vtkRead(mesh_file);
n_nodes = size(mesh.points, 1);
n_cells = size(mesh.cells, 1);
fprintf('Mesh: %d nodes, %d cells\n', n_nodes, n_cells);

%% Load fields
D = load(mat_file);

% Sanity check
assert(size(D.v,1) == n_nodes, ...
    'Node count mismatch: Coarse.vtu=%d, Case_coarse.mat=%d', n_nodes, size(D.v,1));

%% Point data (per node)
mesh = add_point(mesh, 'Fiber',              D.F,           n_nodes);
mesh = add_point(mesh, 'Sheet',              D.F_S,         n_nodes);
mesh = add_point(mesh, 'Normal',             D.F_N,         n_nodes);
mesh = add_point(mesh, 'Transmurality',      D.d3,          n_nodes);
mesh = add_point(mesh, 'Transmurality_RV',   D.Tphi3,       n_nodes);
mesh = add_point(mesh, 'Transmurality_bi',   D.Tphi_bi,     n_nodes);
mesh = add_point(mesh, 'Ventricle',          D.Ventricle,   n_nodes);
mesh = add_point(mesh, 'Epiendo',            D.Epiendo,     n_nodes);
mesh = add_point(mesh, 'Epiendo_RV',         D.Epiendo3,    n_nodes);
mesh = add_point(mesh, 'Apex2Base',          D.a2b,         n_nodes);
mesh = add_point(mesh, 'Apex2Base_cobi',     D.a2b_cobi,    n_nodes);
mesh = add_point(mesh, 'Apex2Base_uvc',      D.a2b_uvc,     n_nodes);
mesh = add_point(mesh, 'Circ_cobi',          D.r,           n_nodes);
mesh = add_point(mesh, 'LVvsRV_cobi',        D.lvrv_cobi,   n_nodes);
mesh = add_point(mesh, 'RV2LV_geo',          D.r2l_geo,     n_nodes);
mesh = add_point(mesh, 'Ant2Post',           D.a2p,         n_nodes);
mesh = add_point(mesh, 'R2L',                D.r2l,         n_nodes);
mesh = add_point(mesh, 'Transmurality_cobi', D.tm_cobi,     n_nodes);
mesh = add_point(mesh, 'Apex2Base_proj',     D.apex_2_base, n_nodes);
mesh = add_point(mesh, 'AHA',                D.aha,         n_nodes);

%% Cell data (per element)
mesh = add_cell(mesh, 'LVvsRV',    D.lvrv,       n_cells);
mesh = add_cell(mesh, 'Label',     D.label_fine, n_cells);
mesh = add_cell(mesh, 'Plug',      D.Plug_tetra, n_cells);
mesh = add_cell(mesh, 'LabelSet2', D.label_set2, n_cells);
mesh = add_cell(mesh, 'AHA',       D.aha,        n_cells);

%% Write
vtkWrite(mesh, out_file);
fprintf('Written: %s\n', out_file);

%% Helper functions
function s = add_point(s, name, val, n)
    if size(val,1) == n
        s.pointData.(name) = double(val);
    else
        warning('Skipping pointData.%s: size %d ~= %d nodes', name, size(val,1), n);
    end
end

function s = add_cell(s, name, val, n)
    if size(val,1) == n
        s.cellData.(name) = double(val);
    else
        warning('Skipping cellData.%s: size %d ~= %d cells', name, size(val,1), n);
    end
end
