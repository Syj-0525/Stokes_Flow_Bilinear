function [K, G] = EssBC_KG(K, G, BV, n_Func)

K_nrows = 3 * n_Func;       % B矩阵行数
for ii = 1: length(BV(:, 1))   % 遍历所有速度边界结点
    NO_u_node = BV(ii, 1);     % 得到速度边界结点编号
    u = BV(ii, 2);
    for jj = 1: K_nrows       % 遍历B矩阵所有行
        G(jj) = G(jj) - K(jj, NO_u_node) * u;
    end
    K(NO_u_node, :) = 0;
    K(:, NO_u_node) = 0;
    K(NO_u_node, NO_u_node) = 1;
    G(NO_u_node) = u;
end

for ii = 1: length(BV(:, 1))   % 遍历所有速度边界结点
    NO_v_node = n_Func + BV(ii, 1);
    v = BV(ii, 3);
    for jj = 1: K_nrows       % 遍历B矩阵所有行
        G(jj) = G(jj) - K(jj, NO_v_node) * v;
    end
    K(NO_v_node, :) = 0;
    K(:, NO_v_node) = 0;
    K(NO_v_node, NO_v_node) = 1;
    G(NO_v_node) = v;
end

%%%%% way2:(强制)速度边界条件实现
% for i = 1:length(BV(:,1))   % 遍历所有速度边界结点
%     NO_u_node = BV(i, 1);
%     u = BV(i, 2);
%     K(NO_u_node, :) = 0;
%     K(NO_u_node, NO_u_node) = 1;
%     G(NO_u_node) = u;
% end
% 
% for i = 1:length(BV(:,1))   % 遍历所有速度边界结点
%     NO_v_node = n_Func + BV(i, 1);
%     v = BV(i, 3);
%     K(NO_v_node, :) = 0;
%     K(NO_v_node, NO_v_node) = 1;
%     G(NO_v_node) = v;
% end
%%%%% way2:(强制)速度边界条件实现

%%% 方腔压力约束点
% K(2 * n_Func + 1, :) = 0;
% K(2 * n_Func + 1, 2 * n_Func + 1) = 1;
% G(2 * n_Func + 1) = 0;
%%% 方腔压力约束点
