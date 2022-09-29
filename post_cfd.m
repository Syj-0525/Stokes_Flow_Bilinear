clc
clear
load solution

%%%%%% 从解向量中获取u, v, p
ux_k_1 = x(1: n_Func);
vy_k_1 = x(1 + n_Func: 2 * n_Func);
p_k_1 = x(1 + 2 * n_Func: 2 * n_Func + n_Func);
%%%%%% 从解向量中获取u, v, p

k = 1; 
for ee = 1: n_el
    for aa = 1: 4
        %%%%%% 生成Tecplot后处理数据
        No_node_post(k, 1) = IEN(ee, aa);
        No_node_post(k, 2) = CtrlPts(IEN(ee, aa), 1);
        No_node_post(k, 3) = CtrlPts(IEN(ee, aa), 2);

        ux_post(k, 1) = ux_k_1(IEN(ee, aa));
        vy_post(k, 1) = vy_k_1(IEN(ee, aa));
        p_post(k, 1) = p_k_1(IEN(ee, aa));  
        %%%%%% 生成Tecplot后处理数据

        k = k + 1;
    end
end

post_data = [No_node_post, ux_post, vy_post, sqrt(ux_post.^2 + vy_post.^2), p_post];
post_data_full = full(post_data);
post_node = [1: 1: size(post_data, 1)];
% post_node2 = reshape(post_node, [n_el/16, 64])';   % 当post_node行太长，对其换行, n_el > 16
%%%%%% 输出Tecplot后处理结果

clear n_el IEN CtrlPts K B N_matrix n_Func NO_node No_node_post p_post      % 清除多余变量
clear viscosity p4 ee aa k p_k_1 ux_k_1 vy_k_1 ux_post vy_post xlabel       % 清除多余变量
save result                                                                 % 存储结果