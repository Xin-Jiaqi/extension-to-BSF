x = linspace(0,1,9); 
y = linspace(0,0.5,5); 

swep_matrix = zeros(length(x)*length(y),2);
count = 1;
for i = 1:length(x)
    for j = 1:length(y)
        swep_matrix(count,:) = ([x(i);y(j)])';
        count = count + 1;
    end
end

% 读取 POSCAR.vesta 文件内容
filename = 'PPdSe2L.vasp';
fileID = fopen(filename, 'r');

% 跳过第一行：标题
fgetl(fileID);

% 读取第二行：缩放系数
scale_factor = str2double(fgetl(fileID));

% 读取第三到第五行：晶格向量
cell_matrix = zeros(3, 3);
for i = 1:3
    cell_matrix(i, :) = sscanf(fgetl(fileID), '%f %f %f');
end

% 读取第六行：元素类型
elements = strtrim(fgetl(fileID));

% 读取第七行：每种元素的原子数量
num_atoms = sscanf(fgetl(fileID), '%d %d %d');

% 读取第八行：坐标类型（Direct）
coordinate_type = strtrim(fgetl(fileID));

% 读取后续行：原子坐标
originalPOSCAR = [];
while ~feof(fileID)
    line = fgetl(fileID);
    if isempty(line)
        break;
    end
    originalPOSCAR = [originalPOSCAR; sscanf(line, '%f %f %f')'];
end

fclose(fileID);

% 构建 POS_inf 文件头信息
POS_inf = sprintf('CONTCAR\n%.1f\n', scale_factor);
for i = 1:3
    POS_inf = [POS_inf, sprintf('%f %f %f\n', cell_matrix(i, :))];
end
POS_inf = [POS_inf, elements, '\n'];
POS_inf = [POS_inf, sprintf('%d %d %d\n', num_atoms), 'Selective dynamics\n', coordinate_type];

i = 1;
for l = 1:length(x)
    for k = 1:length(y)
        newPOSCAR = zeros(size(originalPOSCAR,1),size(originalPOSCAR,2));
        for j = 1:size(originalPOSCAR,1)
            if originalPOSCAR(j,3) > 0.5
                if originalPOSCAR(j,1) + swep_matrix(i,1) >= 1
                    newPOSCAR(j,1) = originalPOSCAR(j,1) + swep_matrix(i,1) - 1;
                else
                    newPOSCAR(j,1) = originalPOSCAR(j,1) + swep_matrix(i,1);
                end
                if originalPOSCAR(j,2) + swep_matrix(i,2) >= 1
                    newPOSCAR(j,2) = originalPOSCAR(j,2) + swep_matrix(i,2) - 1;
                else
                    newPOSCAR(j,2) = originalPOSCAR(j,2) + swep_matrix(i,2);
                end
                newPOSCAR(j,3) = originalPOSCAR(j,3);
            else
                newPOSCAR(j,:) = originalPOSCAR(j,:);
            end
        end
        fileID = fopen([num2str(l),num2str(k),'.vasp'],'w');
        fprintf(fileID,POS_inf);
        for j = 1:size(newPOSCAR,1)
            if newPOSCAR(j,3) > 0.5
                fprintf(fileID,['\n',num2str(newPOSCAR(j,:)),'  F F T']);
            else
                fprintf(fileID,['\n',num2str(newPOSCAR(j,:)),'  F F F']);
            end
        end
        i = i + 1;
    end
end
fclose('all');
