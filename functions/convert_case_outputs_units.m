function convert_case_outputs_units(case_folder, scale, unit_name)
    % Rescale the geometric outputs of one case after the pipeline has run.
    % The pipeline works internally in cm; scale=10 converts outputs to mm.
    %
    % Converted in place:
    %   Coarse.vtu, fine.vtu                         (vtkRead/vtkWrite)
    %   *.vtk legacy ASCII files in the case folder and in EM_Control_*
    %   (labels_final, labelsR, labelsR2, labelspoints, Electrodes,
    %   Activation_points_*)                         (POINTS block only)
    %   ensi*/Fields.ensi.geo coordinates
    %   ensi*/Fields.ensi.Rv_2_lv_geodesic          (geodesic distance)
    %   ensi*/Case*.mat: v and r2l_geo
    %   ensi*/CSVFiles/*_xyz.csv and *_nodefield_cobiveco-rvlv.csv
    % Not converted (solver-specific units): MonoAlg3D/ and Hexa_*.vtk.
    % Input-unit files (input/, labels0.vtk, Closed_final_*.ply) are kept.
    %
    % A marker file prevents a second conversion of the same case.

    if nargin < 3
        unit_name = 'mm';
    end
    marker = fullfile(case_folder, ['units_', unit_name, '.txt']);
    if exist(marker, 'file')
        fprintf('Outputs of %s are already in %s, skipping conversion\n', case_folder, unit_name);
        return
    end

    %volumetric meshes
    for name = {'Coarse.vtu', 'fine.vtu'}
        file = fullfile(case_folder, name{1});
        if exist(file, 'file')
            mesh = vtkRead(file);
            mesh.points = mesh.points .* scale;
            vtkWrite(mesh, file);
        end
    end

    %legacy VTK surfaces and point sets
    legacy = [dir(fullfile(case_folder, '*.vtk')); dir(fullfile(case_folder, 'EM_Control_*', '*.vtk'))];
    for k = 1:numel(legacy)
        if startsWith(legacy(k).name, 'Hexa_') || strcmp(legacy(k).name, 'labels0.vtk')
            continue
        end
        scale_legacy_vtk_points(fullfile(legacy(k).folder, legacy(k).name), scale);
    end

    %field folders
    ensi_dirs = dir(fullfile(case_folder, 'ensi*'));
    ensi_dirs = ensi_dirs([ensi_dirs.isdir]);
    for k = 1:numel(ensi_dirs)
        ensi_path = fullfile(ensi_dirs(k).folder, ensi_dirs(k).name);
        geo = fullfile(ensi_path, 'Fields.ensi.geo');
        if exist(geo, 'file')
            scale_ensi_geo(geo, scale);
        end
        geodesic = fullfile(ensi_path, 'Fields.ensi.Rv_2_lv_geodesic');
        if exist(geodesic, 'file')
            scale_ensi_field(geodesic, scale);
        end

        mats = dir(fullfile(ensi_path, 'Case*.mat'));
        for m = 1:numel(mats)
            file = fullfile(mats(m).folder, mats(m).name);
            data = load(file);
            data.v = data.v .* scale;
            if isfield(data, 'r2l_geo')
                data.r2l_geo = data.r2l_geo .* scale;
            end
            save(file, '-struct', 'data');
        end

        csvs = [dir(fullfile(ensi_path, 'CSVFiles', '*_xyz.csv')); ...
                dir(fullfile(ensi_path, 'CSVFiles', '*_nodefield_cobiveco-rvlv.csv'))];
        for c = 1:numel(csvs)
            file = fullfile(csvs(c).folder, csvs(c).name);
            writematrix(readmatrix(file) .* scale, file);
        end
    end

    fid = fopen(marker, 'w');
    fprintf(fid, 'Geometric outputs scaled by %g to %s (MonoAlg3D/ and Hexa_*.vtk unchanged).\n', scale, unit_name);
    fclose(fid);
    fprintf('Outputs of %s converted to %s\n', case_folder, unit_name);

    %functions
    function scale_legacy_vtk_points(file, s)
        %scale the POINTS block of an ASCII legacy VTK file, keep the rest
        lines = readlines(file);
        if numel(lines) < 3 || ~startsWith(strtrim(lines(3)), 'ASCII')
            warning('Skipping %s: not an ASCII legacy VTK file', file);
            return
        end
        start = find(startsWith(lines, 'POINTS'), 1);
        if isempty(start)
            return
        end
        header = split(strtrim(lines(start)));
        n_values = 3 * str2double(header(2));
        values = zeros(0, 1);
        last = start;
        remainder = strings(0, 1);
        while numel(values) < n_values
            last = last + 1;
            %write_vtk_surf puts the next keyword on the last points line
            [line_values, ~, ~, next] = sscanf(lines(last), '%f');
            values = [values; line_values]; %#ok<AGROW>
            rest = strtrim(extractAfter(lines(last), next - 1));
            if strlength(rest) > 0
                remainder = rest;
            end
        end
        values = reshape(values .* s, 3, []);
        new_points = compose('%.8g %.8g %.8g', values');
        lines = [lines(1:start); new_points; remainder; lines(last+1:end)];
        writelines(lines, file);
    end

    function scale_ensi_geo(file, s)
        %layout written by save_ensi_geo: 6 strings, part id, 2 strings,
        %node count, node ids, then x, y and z arrays (float32)
        fid = fopen(file, 'r+', 'l');
        fseek(fid, 6*80 + 4 + 2*80, 'bof');
        n_nodes = fread(fid, 1, 'uint32');
        offset = ftell(fid) + 4*n_nodes;
        fseek(fid, offset, 'bof');
        xyz = fread(fid, 3*n_nodes, 'float32');
        fseek(fid, offset, 'bof');
        fwrite(fid, xyz .* s, 'float32');
        fclose(fid);
    end

    function scale_ensi_field(file, s)
        %layout written by save_ensi_field: 2 strings, part id, 1 string,
        %then the node values (float32)
        fid = fopen(file, 'r+', 'l');
        offset = 2*80 + 4 + 80;
        fseek(fid, offset, 'bof');
        values = fread(fid, Inf, 'float32');
        fseek(fid, offset, 'bof');
        fwrite(fid, values .* s, 'float32');
        fclose(fid);
    end

end
