function convert_poscar_cartesian_to_direct(input_file, output_file)
    % 读取POSCAR文件
    fileID = fopen(input_file, 'r');
    
    % 读取文件头
    title = fgetl(fileID); % 第一行：标题信息
    scale_factor = str2double(fgetl(fileID)); % 第二行：缩放因子
    lattice_vectors = zeros(3, 3); % 初始化晶格矩阵
    
    % 读取晶格参数
    for i = 1:3
        lattice_vectors(i, :) = fscanf(fileID, '%f %f %f', [1, 3]);
    end
    
    % 读取原子种类和数量
    fgetl(fileID);
    atom_types_line = fgetl(fileID); % 第6行：原子种类
    atom_counts_line = fgetl(fileID); % 第7行：每种原子数量
    atom_types = strsplit(strtrim(atom_types_line)); % 提取原子种类
    atom_counts = str2num(atom_counts_line); %#ok<ST2NM> % 提取原子数量
    total_atoms = sum(atom_counts); % 总原子数
    
    % 读取坐标类型
    coord_type = strtrim(fgetl(fileID)); % 第8行：坐标类型
    
    % 检查坐标是否是笛卡尔形式
    if ~strcmpi(coord_type, 'Cartesian')
        error('The input POSCAR file is not in Cartesian coordinate format.');
    end
    
    % 读取原子坐标
    cartesian_coords = fscanf(fileID, '%f %f %f', [3, total_atoms])'; % 读取坐标
    
    % 关闭文件
    fclose(fileID);
    
    % 计算分数坐标
    fractional_coords = cartesian_coords / lattice_vectors / scale_factor;
    
    % 写入新的POSCAR文件
    fileID = fopen(output_file, 'w');
    fprintf(fileID, '%s\n', title);
    fprintf(fileID, '1.0\n'); % 缩放因子始终设置为1.0，保持原始晶格不变
    
    for i = 1:3
        fprintf(fileID, '%20.12f %20.12f %20.12f\n', lattice_vectors(i, :));
    end
    
    fprintf(fileID, '%s\n', atom_types_line); % 原子种类
    fprintf(fileID, '%s\n', atom_counts_line); % 原子数量
    fprintf(fileID, 'Direct\n'); % 转换后坐标类型
    
    % 写入分数坐标
    for i = 1:total_atoms
        fprintf(fileID, '%20.12f %20.12f %20.12f\n', fractional_coords(i, :));
    end
    
    fclose(fileID);
    
    fprintf('Conversion completed. The Direct coordinates have been saved to %s.\n', output_file);
end