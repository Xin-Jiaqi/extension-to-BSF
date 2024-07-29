function create_bilayer_vasp_files(filename, matfile)
    % 加载.mat文件中的变换矩阵
    loaded_data = load(matfile);

    % 获取.mat文件中的所有字段名
    matrix_fields = fieldnames(loaded_data);

    % 初始化变换矩阵的 cell 数组
    transforms = cell(length(matrix_fields), 1);

    % 提取每个变换矩阵并存储到 cell 数组中
    for idx = 1:length(matrix_fields)
        transforms{idx} = loaded_data.(matrix_fields{idx});
    end

    % 确保文件名没有扩展名
    [~, name, ~] = fileparts(filename);

    % 创建用于存储调整后的目录
    bilayers_dir = 'bilayers';
    if ~exist(bilayers_dir, 'dir')
        mkdir(bilayers_dir);
    end

    % 读取POSCAR文件
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件: %s', filename);
    end

    % 读取前6行，POSCAR文件的前6行包含其他信息
    header = cell(6, 1);
    for i = 1:6
        header{i} = fgetl(fid);
    end

    % 第7行包含原子种类的数量
    atom_counts_line = fgetl(fid);

    % 将第7行的内容转换为数字数组
    atom_counts = sscanf(atom_counts_line, '%d');

    % 计算总原子数
    atom_count = sum(atom_counts);

    % 输出总原子数
    fprintf('POSCAR文件中的总原子数: %d\n', atom_count);

    % 读取第8行，通常为选择坐标系类型的一行
    % 这行通常是 "Direct" 或者 "Cartesian"
    coordinate_type = strtrim(fgetl(fid));

    % 读取从第9行开始的原子坐标，存储到矩阵中
    coordinates = zeros(atom_count, 3);  % 初始化原子坐标矩阵
    for i = 1:atom_count
        coord_line = fgetl(fid);
        coordinates(i, :) = sscanf(coord_line, '%f %f %f');
    end

    % 关闭文件
    fclose(fid);

    % 判断是否需要从 Cartesian 转换为 Direct
    if strcmpi(coordinate_type, 'Cartesian')
        disp('坐标为笛卡尔坐标，将转换为分数坐标 (Direct)');

        % 从文件中读取晶格参数
        lattice_vectors = zeros(3, 3);
        fid = fopen(filename, 'r');
        fgetl(fid); % 跳过第一行标题
        scale_factor = str2double(fgetl(fid)); % 第二行: 缩放因子

        for i = 1:3
            lattice_vectors(i, :) = fscanf(fid, '%f %f %f', [1, 3]);
        end
        fclose(fid);

        % 计算分数坐标
        coordinates = (coordinates / lattice_vectors) / scale_factor;
    elseif ~strcmpi(coordinate_type, 'Direct')
        error('未知的坐标系类型: %s', coordinate_type);
    end

    % 找出原始POSCAR文件中第三列的最小值
    min_z = min(coordinates(:, 3));

    % 计算0.1与最小值的差
    delta_z_orig = 0.1 - min_z;

    % 将差值加到原始坐标的第三列上
    modified_orig_coordinates = coordinates;
    modified_orig_coordinates(:, 3) = modified_orig_coordinates(:, 3) + delta_z_orig;

    % 确保调整后的原始坐标也满足周期性边界条件
    for i = 1:size(modified_orig_coordinates, 1)
        if modified_orig_coordinates(i, 3) < 0
            modified_orig_coordinates(i, 3) = modified_orig_coordinates(i, 3) + 1;
        elseif modified_orig_coordinates(i, 3) > 1
            modified_orig_coordinates(i, 3) = modified_orig_coordinates(i, 3) - 1;
        end
    end

    % 遍历每个变换矩阵并应用到原子坐标
    for idx = 1:length(transforms)
        trans = transforms{idx};  % 当前变换矩阵

        % 应用变换矩阵到原子坐标
        transformed_coordinates = (trans * coordinates')';

        % 对转换后的坐标应用周期性边界条件处理
        for i = 1:size(transformed_coordinates, 1)
            for j = 1:size(transformed_coordinates, 2)
                if transformed_coordinates(i, j) < 0
                    transformed_coordinates(i, j) = transformed_coordinates(i, j) + 1;
                elseif transformed_coordinates(i, j) > 1
                    transformed_coordinates(i, j) = transformed_coordinates(i, j) - 1;
                end
            end
        end

        % 获取当前变换矩阵的名称
        matrix_name = matrix_fields{idx};

        %% 处理 bilayers 的额外步骤 %%
        % 找出第三列的最大值
        max_z = max(transformed_coordinates(:, 3));

        % 计算0.9与最大值的差
        delta_z = 0.9 - max_z;

        % 将差值加到第三列上
        bilayer_coordinates = transformed_coordinates;
        bilayer_coordinates(:, 3) = bilayer_coordinates(:, 3) + delta_z;

        % 确保调整后的坐标也满足周期性边界条件
        for i = 1:size(bilayer_coordinates, 1)
            if bilayer_coordinates(i, 3) < 0
                bilayer_coordinates(i, 3) = bilayer_coordinates(i, 3) + 1;
            elseif bilayer_coordinates(i, 3) > 1
                bilayer_coordinates(i, 3) = bilayer_coordinates(i, 3) - 1;
            end
        end

        % 构造 bilayer 输出文件名
        bilayer_output_filename = fullfile(bilayers_dir, sprintf('%s_%s_bilayer.vasp', name, matrix_name));

        % 打开新的 bilayer 文件用于写入
        fid_bilayer_out = fopen(bilayer_output_filename, 'w');
        if fid_bilayer_out == -1
            error('无法创建文件: %s', bilayer_output_filename);
        end

        % 写入POSCAR文件的前6行
        for i = 1:5
            fprintf(fid_bilayer_out, '%s\n', header{i});
        end

        % 将原始POSCAR的第6行和第7行直接写在新 bilayer 文件对应行的后面
        fprintf(fid_bilayer_out, '%s %s\n', header{6}, header{6}); % 第6行
        fprintf(fid_bilayer_out, '%s %s\n', atom_counts_line, atom_counts_line); % 第7行

        % 写入坐标系类型
        fprintf(fid_bilayer_out, '%s\n', 'Direct');

        % 写入调整后的坐标
        for i = 1:atom_count
            fprintf(fid_bilayer_out, '%12.8f %12.8f %12.8f\n', bilayer_coordinates(i, :));
        end

        % 写入处理过的原始POSCAR坐标
        for i = 1:atom_count
            fprintf(fid_bilayer_out, '%12.8f %12.8f %12.8f\n', modified_orig_coordinates(i, :));
        end

        % 关闭 bilayer 输出文件
        fclose(fid_bilayer_out);

        % 输出 bilayer 结果文件名
        fprintf('生成的 bilayer 文件: %s\n', bilayer_output_filename);
    end
end

% 调用函数并传入POSCAR文件名和mat文件名
% 例如，create_bilayer_vasp_files('POSCAR', 'transform_matrices.mat');