clc
clear
clf
%%%%%%% 区域几何尺寸及网格划分参数
H = 1; % 区域总高（m）
L = 2; % 区域总长（m）
nel_x = 8; % 水平方向单元数量
nel_y = 8; % 竖直方向单元数量
%%%%%%% 区域几何尺寸及网格划分参数

%%%%%%% 总单元及总结点数
n_el = nel_x * nel_y; % 总单元数
n_Func = (nel_x + 1) * (nel_y + 1); % 总结点数
%%%%%%% 总单元及总结点数

%%%%%%% 单元间距
dx = L/nel_x; % 水平方向网格间距
dy = H/nel_y; % 数值方向网格间距
%%%%%%% 单元间距

%%%%%%% 节点分布编号
NO_node = zeros(nel_x + 1, nel_y + 1);
for i = 1:1:nel_y + 1
    for j = 1:1:nel_x + 1
    NO_node(i, j) = (i - 1) * (nel_x + 1) + j;
    end
end
%%%%%%% 节点分布编号

%%%%%%% CtrlPts(双线性单元结点坐标数据)
%%%%%%% 四边形双线性单元CtrlPts生成
CtrlPts = zeros(n_Func, 2);
for i = 1:1:nel_y + 1
    for j = 1:1:nel_x + 1
        CtrlPts(NO_node(i, j), 1) = dx * (j - 1);
        CtrlPts(NO_node(i, j), 2) = dy * (i - 1);
    end
end
%%%%%%% 四边形双线性单元CtrlPts生成

%%%%%%% IEN(四边形双线性单元IEN)
%%%%%%% 四边形双线性单元IEN生成
IEN = zeros(n_el, 4);
ee = 0;
for i = 1:nel_y
    for j = 1:nel_x
        ee = ee + 1;
        IEN(ee,:) = [NO_node(i, j), NO_node(i, j + 1), NO_node(i + 1, j + 1), NO_node(i + 1, j),];
    end
end
%%%%%%% 四边形双线性单元IEN生成

%%%%%%% NBC(边界节点数据)：边界节点编号
%%%%%%% NBC数据生成
NBC1 = NO_node(1,:);
NBC2 = NO_node(:, nel_x + 1);
NBC3 = NO_node(nel_y + 1, :);
NBC4 = NO_node(:,1);
%%%%%%% NBC数据生成

%%%%%%% EBC边界单元数据：边界单元编号，边界单元上处于边界的边的编号，nx, ny
%%%%%%% EBC数据生成
thetax1 = pi/2;      % 1号边界外法线方向与x轴夹角
thetay1 = pi;        % 1号边界外法线方向与y轴夹角
thetax2 = 0;         % 2号边界外法线方向与x轴夹角
thetay2 = pi/2;      % 2号边界外法线方向与y轴夹角
thetax3 = pi/2;      % 3号边界外法线方向与x轴夹角
thetay3 = 0;         % 3号边界外法线方向与y轴夹角
thetax4 = pi;        % 4号边界外法线方向与x轴夹角
thetay4 = pi/2;      % 4号边界外法线方向与y轴夹角
EBC1 = [ (1: nel_x)', ones(size([1: nel_x]')),  ones(nel_x, 1) * cos(thetax1), ones(nel_x, 1) * cos(thetay1)];
EBC2 = [ (nel_x: nel_x: nel_y * nel_x)', 2 * ones(size([1:nel_y]')), ones(nel_y, 1) * cos(thetax2), ones(nel_y, 1) * cos(thetay2)];
EBC3 = [ (nel_x * (nel_y - 1) + 1: nel_x * nel_y)', 3 * ones(size([1: nel_x]')), ones(nel_x, 1) * cos(thetax3), ones(nel_x, 1) * cos(thetay3)];
EBC4 = [ (1:nel_x:nel_x * (nel_y - 1) + 1)', 4 * ones(size([1:nel_y]')), ones(nel_y, 1) * cos(thetax4), ones(nel_y, 1) * cos(thetay4)];
%%%%%%% EBC数据生成

%%%%%%% 调用四边形网格绘制程序
% plot_mesh(IEN, CtrlPts);
%%%%%%% 调用四边形网格绘制程序

%%%%%%% 清楚多余变量
clear dx dy H L nel_x nel_y i j k ee
clear theta No_node
clear thetax1 thetax2 thetax3 thetax4
clear thetay1 thetay2 thetay3 thetay4
%%%%%%% 清楚多余变量

%%%%%%% 存储网格数据
save msh
% save ('../source/matlab_repos/Computational mechanics/Project_Final/msh')
%%%%%%% 存储网格数据