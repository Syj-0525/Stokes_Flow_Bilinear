clc
clear
format long
%%%%%% 读取网格数据
load msh
%%%%%% 读取网格数据

%%%%%% 设定流体黏度
viscosity = 1000;
%%%%%% 设定流体黏度

%%%%%% 设定边界条件
u1 = 0; v1 = 0;     % 设定1边界速度
u3 = 0.01; v3 = 0;  % 设定3边界速度

%%% Case1 槽道流
%%%%% 设定边界条件数据BP与BV
%%% BV：边界节点编号，该节点u, 该节点v
BV1 =[ NBC1', u1 * ones(size(NBC1))', v1 * ones(size(NBC1))' ];
BV3 =[ NBC3', u3 * ones(size(NBC3))', v3 * ones(size(NBC3))' ];
BV = [ BV1; BV3 ];
P2 = 0;             % 设定边界压强
P4 = 1000;          % 设定边界压强
%%% BP：边界单元编号，边界单元上处于边界的边的编号，nx, ny，边界单元上处于边界的边上第一节点压力 第二节点压力
BP2 = [ EBC2 ,ones(size(EBC2(:,1))) * P2, ones(size(EBC2(:,1))) * P2 ];
BP4 = [ EBC4 ,ones(size(EBC4(:,1))) * P4, ones(size(EBC4(:,1))) * P4 ];
BP = [ BP2; BP4 ];
%%%%%% 设定边界条件数据BP与BV

%%% Case2
% %%% 方腔边界
% BV1 =[ NBC1', u1 * ones(size(NBC1))', v1 * ones(size(NBC1))'];
% BV3 =[ NBC3', u3 * ones(size(NBC3))', v3 * ones(size(NBC3))'];
% u2 = 0; v2 = 0;     % 设定2边界速度
% u4 = 0; v4 = 0;     % 设定4边界速度
% BV2 =[ NBC2, u2 * ones(size(NBC2)), v2 * ones(size(NBC2))];
% BV4 =[ NBC4, u4 * ones(size(NBC4)), v4 * ones(size(NBC4))];
% BV = [BV1; BV3; BV2; BV4];
% P2 = 0;             % 设定方腔边界压强
% P4 = 0;             % 设定方腔边界压强
% BP2 = [ EBC2 ,ones(size(EBC2(:,1))) * P2, ones(size(EBC2(:,1))) * P2 ];
% BP4 = [ EBC4 ,ones(size(EBC4(:,1))) * P4, ones(size(EBC4(:,1))) * P4 ];
% BP = [BP2; BP4];
% %%% 方腔边界

clear BV1 BV3 BV2 BV4 BP1 BP2 BP3 BP4 P2 P4
clear EBC1 EBC2 EBC3 EBC4
clear NBC1 NBC2 NBC3 NBC4
clear u1 v1 u3 v3 u2 v2 u4 v4
%%%%%% 设定边界条件

%%%%%% 获取全局刚度矩阵
[K, G] = GAssem_KG(n_Func, n_el, IEN, CtrlPts, viscosity, BP, BV);
%%%%%% 获取全局刚度矩阵

%%%%%% 求解线性方程组
x = Linear_solver(K, G);
%%%%%% 求解线性方程组

%%%%%% 清理变量
clear BP BV viscosity
%%%%%% 清理变量

% 存储解向量
save solution
% 存储解向量



