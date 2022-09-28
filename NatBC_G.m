function[G] = NatBC_G(BP, CtrlPts, IEN, n_Func)

% F1 = zeros(n_Func, 1);
% F2 = zeros(n_Func, 1);
% F3 = zeros(n_Func, 1);

F1 = sparse(n_Func, 1);
F2 = sparse(n_Func, 1);
F3 = sparse(n_Func, 1);


eleCtrlPts = zeros(4 , 2);  % 4：一个单元4个节点，2：x y
%%%%%% 计算各单元荷载(子块)矩阵并装配为整体荷载的子块矩阵
for ii = 1:length(BP(:, 1))   % 遍历所有压强边界单元
    for aa = 1:4
        eleCtrlPts(aa,:) = CtrlPts(IEN(BP(ii, 1), aa), :);   % 第i个压强边界单元的速度结点坐标数据
    end
    [Fe1, Fe2] = LocAssem_Fe(eleCtrlPts, BP(ii, :));   % 组装 Fe

    %%%%%% 装配F1, F2
    for aa = 1:4
        F1(IEN(BP(ii, 1), aa), 1) = F1(IEN(BP(ii, 1), aa), 1) + Fe1(aa, 1);
        F2(IEN(BP(ii, 1), aa), 1) = F2(IEN(BP(ii, 1), aa), 1) + Fe2(aa, 1);
    end
    %%%%%% 装配F1, F2
end
%%%%%% 计算各单元荷载(子块)矩阵并装配为整体荷载的子块矩阵

G = [ -F1; -F2; -F3 ];