% Debug script: compare labelling axes between UKBB and NUH
% Run from InSilicoHeartGen root directory

clear
addpath(genpath(fullfile(pwd, 'functions')));
addpath(genpath(fullfile(pwd, 'dependencies')));

%% Load NUH surface
nuh_path = fullfile(pwd, 'outputs', 'NUH_test', '1');
cd(nuh_path);
sur_nuh = vtkRead('labels_final.vtk');

%% Load UKBB surface (if available)
ukbb_path = fullfile(pwd, '..', '..', 'UKBB', '1');
if exist(fullfile(ukbb_path, 'labels_final.vtk'), 'file')
    cd(ukbb_path);
    sur_ukbb = vtkRead('labels_final.vtk');
    has_ukbb = true;
else
    has_ukbb = false;
    fprintf('No UKBB labels_final.vtk found at %s\n', ukbb_path);
end

%% Analyze NUH
fprintf('\n=== NUH Analysis ===\n');
analyze_labels(sur_nuh, 'NUH');

if has_ukbb
    fprintf('\n=== UKBB Analysis ===\n');
    analyze_labels(sur_ukbb, 'UKBB');
end

function analyze_labels(sur, name)
    node_surf = sur.points;
    face_surf = sur.cells;
    labelfinal = sur.cellData;

    % Get the cell data field name
    fnames = fieldnames(labelfinal);
    labelfinal = labelfinal.(fnames{1});

    centroid = meshcentroid(node_surf, face_surf);

    % Print label distribution
    unique_labels = unique(labelfinal);
    fprintf('Labels present: ');
    for i = 1:length(unique_labels)
        fprintf('%d(%d) ', unique_labels(i), sum(labelfinal == unique_labels(i)));
    end
    fprintf('\n');

    % Label meanings:
    % 1=RV_epi, 2=LV_endo, 3=RV_endo, 4=LV_epi
    % 6=septum, 9=PV_endo, 10=MV, 12=LV_apex, 13=AV
    % 14=TV_endo, 18=RV_apex, 19=PV_epi, 24=TV_epi

    % Check PV location (label 9)
    if any(labelfinal == 9)
        pv_centroid = mean(centroid(labelfinal == 9, :));
        fprintf('PV center (label 9): [%.2f, %.2f, %.2f], %d cells\n', pv_centroid, sum(labelfinal == 9));
    else
        fprintf('NO PV (label 9) found!\n');
    end

    % Check TV location (label 14)
    if any(labelfinal == 14)
        tv_centroid = mean(centroid(labelfinal == 14, :));
        fprintf('TV center (label 14): [%.2f, %.2f, %.2f], %d cells\n', tv_centroid, sum(labelfinal == 14));
    else
        fprintf('NO TV (label 14) found!\n');
    end

    % Check RV apex (label 18)
    if any(labelfinal == 18)
        apex_rv = centroid(labelfinal == 18, :);
        fprintf('RV apex (label 18): [%.2f, %.2f, %.2f]\n', apex_rv);
    end

    % Check LV apex (label 12)
    if any(labelfinal == 12)
        apex_lv = centroid(labelfinal == 12, :);
        fprintf('LV apex (label 12): [%.2f, %.2f, %.2f]\n', apex_lv);
    end

    % AV location (label 13)
    if any(labelfinal == 13)
        av_centroid = mean(centroid(labelfinal == 13, :));
        fprintf('AV center (label 13): [%.2f, %.2f, %.2f], %d cells\n', av_centroid, sum(labelfinal == 13));
    end

    % MV location (label 10)
    if any(labelfinal == 10)
        mv_centroid = mean(centroid(labelfinal == 10, :));
        fprintf('MV center (label 10): [%.2f, %.2f, %.2f], %d cells\n', mv_centroid, sum(labelfinal == 10));
    end

    % Compute key axes (same as Ventricular_Labelling.m)
    TR = triangulation(double(face_surf), double(node_surf));
    face_normals = faceNormal(TR);

    % RV axis (from AV normals, label 13)
    if any(labelfinal == 13)
        RVaxis = median(face_normals(labelfinal == 13, :));
        fprintf('RVaxis (median AV normals): [%.4f, %.4f, %.4f]\n', RVaxis);

        % computeLongAxis for RV
        surfRV.cells = face_surf(labelfinal == 3 | labelfinal == 9 | labelfinal == 14 | labelfinal == 18, :);
        % Use original label 3 (RV endo) cells
        surfRV.cells = face_surf(any(labelfinal == [3], 2), :);
        surfRV.points = node_surf;
        surfRV.cellTypes = uint8(ones(size(surfRV.cells, 1), 1) * 5);

        longAxRV = computeLongAxis(surfRV, RVaxis);
        fprintf('longAxRV (computeLongAxis): [%.4f, %.4f, %.4f]\n', longAxRV);
    end

    % LV2RV direction
    lv_center = mean(centroid(labelfinal == 2, :));
    rv_center = mean(centroid(labelfinal == 3, :));
    LV2RV = rv_center - lv_center;
    fprintf('LV2RV (RV_center - LV_center): [%.4f, %.4f, %.4f]\n', LV2RV);

    % Cross product for anterior-posterior
    RVaxispulm_cross = cross(LV2RV, RVaxis);
    RVaxispulm_cross_norm = RVaxispulm_cross / norm(RVaxispulm_cross);
    fprintf('cross(LV2RV, RVaxis) normalized: [%.4f, %.4f, %.4f]\n', RVaxispulm_cross_norm);

    % computeLongAxis on the cross product
    RVaxispulm_opt = computeLongAxis(surfRV, RVaxispulm_cross);
    fprintf('computeLongAxis(surfRV, cross): [%.4f, %.4f, %.4f]\n', RVaxispulm_opt);

    % Angle between cross product and optimized result
    angle_deg = acosd(abs(dot(RVaxispulm_cross_norm, RVaxispulm_opt)));
    fprintf('Angle between cross and optimized: %.1f degrees\n', angle_deg);

    % Where does max heightRV2 point to?
    heightRV2_cross = double(node_surf * RVaxispulm_cross_norm');
    heightRV2_opt = double(node_surf * RVaxispulm_opt');

    [~, max_cross_id] = max(heightRV2_cross);
    [~, max_opt_id] = max(heightRV2_opt);

    fprintf('Max heightRV2 point (cross): [%.2f, %.2f, %.2f]\n', node_surf(max_cross_id, :));
    fprintf('Max heightRV2 point (optimized): [%.2f, %.2f, %.2f]\n', node_surf(max_opt_id, :));
end
