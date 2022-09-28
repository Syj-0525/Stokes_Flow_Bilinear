function [K, G] = GAssem_KG(n_Func, n_el, IEN, CtrlPts, viscosity, BP, BV)

%%%%%% 初始化总体方程的各个子块矩阵
% B1 = zeros(n_Func, n_Func);
% B2 = zeros(n_Func, n_Func);
% D11 = zeros(n_Func, n_Func);
% D12 = zeros(n_Func, n_Func);
% D21 = zeros(n_Func, n_Func);
% D22 = zeros(n_Func, n_Func);
% C1 = zeros(n_Func, n_Func);
% C2 = zeros(n_Func, n_Func);
% D = zeros(n_Func, n_Func)

C1 = sparse(n_Func, n_Func);   %采用稀疏阵（稀疏阵赋值慢）
C2 = sparse(n_Func, n_Func);
A11 = sparse(n_Func, n_Func);
A12 = sparse(n_Func, n_Func);
A21 = sparse(n_Func, n_Func);
A22 = sparse(n_Func, n_Func);
B1 = sparse(n_Func, n_Func);
B2 = sparse(n_Func, n_Func);
D = sparse(n_Func, n_Func);
%%%%%% 初始化总体方程的各个子块矩阵

eleCtrlPts = zeros(4 , 2);  % 4：一个单元4个节点，2：x y
%%%%%% 计算各单元刚度(子块)矩阵并装配为整体刚度的子块矩阵
for ii = 1: n_el   % 遍历所有单元
    for aa = 1: 4  % 遍历单元结点
        eleCtrlPts(aa,:) = CtrlPts(IEN(ii, aa), :);   % 第i个速度单元的结点坐标数据
    end
    [Ce1, Ce2] = LocAssem_Ce(eleCtrlPts);   % 组装 Ce
    [Ae11, Ae12, Ae21, Ae22] = LocAssem_Ae(eleCtrlPts, viscosity);   % 组装 Ae
    [Be1, Be2] = LocAssem_Be(eleCtrlPts);   % 组装 Be

    %%%%%% 装配B1, B2, A11, A12, A21, A22, C1, C2
    for aa = 1:4
        for bb = 1:4
            C1(IEN(ii, aa), IEN(ii, bb)) = C1(IEN(ii, aa), IEN(ii, bb)) + Ce1(aa, bb);
            C2(IEN(ii, aa), IEN(ii, bb)) = C2(IEN(ii, aa), IEN(ii, bb)) + Ce2(aa, bb);
            A11(IEN(ii, aa), IEN(ii, bb)) = A11(IEN(ii, aa), IEN(ii, bb)) + Ae11(aa, bb);
            A12(IEN(ii, aa), IEN(ii, bb)) = A12(IEN(ii, aa), IEN(ii, bb)) + Ae12(aa, bb);
            A21(IEN(ii, aa), IEN(ii, bb)) = A21(IEN(ii, aa), IEN(ii, bb)) + Ae21(aa, bb);
            A22(IEN(ii, aa), IEN(ii, bb)) = A22(IEN(ii, aa), IEN(ii, bb)) + Ae22(aa, bb);  
            B1(IEN(ii, aa), IEN(ii, bb)) = B1(IEN(ii, aa), IEN(ii, bb)) + Be1(aa, bb);
            B2(IEN(ii, aa), IEN(ii, bb)) = B2(IEN(ii, aa), IEN(ii, bb)) + Be2(aa, bb);  
        end
    end
    %%%%%% 装配B1, B2, A11, A12, A21, A22, C1, C2
end
%%%%%% 计算各单元刚度(子块)矩阵并装配为整体刚度的子块矩阵

%%%%%% 组合为整体矩阵方程
K =  [A11 A12 -B1
      A21 A22 -B2
      C1  C2   D];   
%%%%%% 组合为整体矩阵方程

%%%%%% 自然边界条件实现
[G] = NatBC_G(BP, CtrlPts, IEN, n_Func);
%%%%%% 自然边界条件实现

%%%%%% (强制)速度边界条件实现
[K, G] = EssBC_KG(K, G, BV, n_Func);
%%%%%% (强制)速度边界条件实现


