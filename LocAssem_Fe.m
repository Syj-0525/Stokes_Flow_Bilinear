function [Fe1, Fe2] = LocAssem_Fe(eleCtrlPts, eleBP)
%%%%%%% 初始化Fe1和Fe2
Fe1 = zeros(4, 1);
Fe2 = zeros(4, 1);
%%%%%%% 初始化Fe1和Fe2

%%%%%%% 高斯积分点和权重
nqp = 10;
[gp, gw] = Gauss(nqp, -1, 1);
% gp = [0.9324695142031521, 0.6612093864662645, 0.2386191860831969, -0.9324695142031521, -0.6612093864662645, -0.2386191860831969];
% gw = [0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.1713244923791704, 0.3607615730481386, 0.4679139345726910];
kesi = gp;
ita = gp;
%%%%%%% 高斯积分点和权重

%%%%%%% (压强)自然边界条件数据
No_p_edge = eleBP(1, 2);    % 边界单元的边界编号
nx = eleBP(1, 3);           % 边界边外法向与x轴余弦
ny = eleBP(1, 4);           % 边界边外法向与y轴余弦
p_val1 = eleBP(1, 5);     % 边界单元边第一点处压力值
p_val2 = eleBP(1, 6);     % 边界单元边第二点处压力值
%%%%%%% (压强)自然边界条件数据

%%%%%%% Fe1, Fe2的数值积分
%%%%%%% 当压力施加于单元第1边时，Fe1和Fe2的计算
if No_p_edge == 1
    for ii = 1: nqp
        %%%%%%% ita = -1时，形函数及其对kesi的导数
        ita = -1;
        %%%%%%% 形函数
        R = [1/4 * (1 - kesi(ii)) * (1 - ita)
             1/4 * (1 + kesi(ii)) * (1 - ita)
             1/4 * (1 + kesi(ii)) * (1 + ita)
             1/4 * (1 - kesi(ii)) * (1 + ita)];
        %%%%%%% 形函数
      
        P = [p_val1, p_val2, 0, 0]';

        R_kesi = [ - 1/4 * (1 - ita)
                     1/4 * (1 - ita)
                     1/4 * (1 + ita)
                   - 1/4 * (1 + ita) ];
        %%%%%%% ita = -1时，形函数及其对对kesi的导数

        %%%%%% Jacobi相关计算
        dx_dkesi = R_kesi' * eleCtrlPts(:,1);
        dy_dkesi = R_kesi' * eleCtrlPts(:,2);
        det_Jacobi = sqrt(dx_dkesi^2 + dy_dkesi^2);
        %%%%%% Jacobi相关计算

%         Fe1 = Fe1 + gw(ii) * R' * P * nx * R * det_Jacobi;
%         Fe2 = Fe2 + gw(ii) * R' * P * ny * R * det_Jacobi;
        for aa = 1: 4
            Na = R(aa);
            Fe1(aa, 1) = Fe1(aa, 1) + gw(ii) * Na * (R' * P) * nx * det_Jacobi;
            Fe2(aa, 1) = Fe2(aa, 1) + gw(ii) * Na * (R' * P) * ny * det_Jacobi;
        end
    end
end
%%%%%%% 当压力施加于单元第1边时，Fe1和Fe2的计算

%%%%%%% 当压力施加于单元第2边时，Fe1和Fe2的计算
if No_p_edge == 2
    for jj = 1: nqp
        %%%%%%% kesi = 1时，形函数及其对ita的导数
        kesi = 1;
        R = [1/4 * (1 - kesi) * (1 - ita(jj))
             1/4 * (1 + kesi) * (1 - ita(jj))
             1/4 * (1 + kesi) * (1 + ita(jj))
             1/4 * (1 - kesi) * (1 + ita(jj))];
        P = [0, p_val1, p_val2, 0]';

        R_ita = [ - 1/4 * (1 - kesi)
                  - 1/4 * (1 + kesi)
                    1/4 * (1 + kesi)
                    1/4 * (1 - kesi) ];
        %%%%%%% kesi = 1时，形函数及其对ita的导数

        %%%%%% Jacobi相关计算
        dx_dita = R_ita' * eleCtrlPts(:,1);
        dy_dita = R_ita' * eleCtrlPts(:,2);
        det_Jacobi = sqrt(dx_dita^2 + dy_dita^2);
        %%%%%% Jacobi相关计算

%         Fe1 = Fe1 + gw(jj) * R' * P * nx * R * det_Jacobi;
%         Fe2 = Fe2 + gw(jj) * R' * P * ny * R * det_Jacobi;
        for aa = 1: 4
            Na = R(aa);
            Fe1(aa, 1) = Fe1(aa, 1) + gw(jj) * Na * (R' * P) * nx * det_Jacobi;
            Fe2(aa, 1) = Fe2(aa, 1) + gw(jj) * Na * (R' * P) * ny * det_Jacobi;
        end
    end
end
%%%%%%% 当压力施加于单元第2边时，Fe1和Fe2的计算

%%%%%%% 当压力施加于单元第3边时，Fe1和Fe2的计算
if No_p_edge == 3
    for ii = 1: nqp
        %%%%%%% ita = 1时，形函数及其对kesi的导数
        ita = 1;
        R = [1/4 * (1 - kesi(ii)) * (1 - ita)
             1/4 * (1 + kesi(ii)) * (1 - ita)
             1/4 * (1 + kesi(ii)) * (1 + ita)
             1/4 * (1 - kesi(ii)) * (1 + ita)];
        P = [0, 0, p_val1, p_val2]';

        R_kesi = [ - 1/4 * (1 - ita)
                     1/4 * (1 - ita)
                     1/4 * (1 + ita)
                   - 1/4 * (1 + ita) ];
        %%%%%%% ita = 1时，形函数及其对kesi的导数

        %%%%%% Jacobi相关计算
        dx_dkesi = R_kesi' * eleCtrlPts(:,1);
        dy_dkesi = R_kesi' * eleCtrlPts(:,2);
        det_Jacobi = sqrt(dx_dkesi^2 + dy_dkesi^2);
        %%%%%% Jacobi相关计算

%         Fe1 = Fe1 + gw(ii) * R' * P * nx * R * det_Jacobi;
%         Fe2 = Fe2 + gw(ii) * R' * P * ny * R * det_Jacobi;
        for aa = 1: 4
            Na = R(aa);
            Fe1(aa, 1) = Fe1(aa, 1) + gw(ii) * Na * (R' * P) * nx * det_Jacobi;
            Fe2(aa, 1) = Fe2(aa, 1) + gw(ii) * Na * (R' * P) * ny * det_Jacobi;
        end
    end
end
%%%%%%% 当压力施加于单元第3边时，Fe1和Fe2的计算

%%%%%%% 当压力施加于单元第4边时，Fe1和Fe2的计算
if No_p_edge == 4
    for jj = 1: nqp
        %%%%%%% kesi = -1时，形函数及其对kesi的导数
        kesi = -1;
        R = [1/4 * (1 - kesi) * (1 - ita(jj))
             1/4 * (1 + kesi) * (1 - ita(jj))
             1/4 * (1 + kesi) * (1 + ita(jj))
             1/4 * (1 - kesi) * (1 + ita(jj))];
        P = [p_val1, 0, 0, p_val2]';

        R_ita = [ - 1/4 * (1 - kesi)
                  - 1/4 * (1 + kesi)
                    1/4 * (1 + kesi)
                    1/4 * (1 - kesi) ];
        %%%%%%% kesi = -1时，形函数及其对kesi的导数

        %%%%%% Jacobi相关计算
        dx_dita = R_ita' * eleCtrlPts(:,1);
        dy_dita = R_ita' * eleCtrlPts(:,2);
        det_Jacobi = sqrt(dx_dita^2 + dy_dita^2);
        %%%%%% Jacobi相关计算

%         Fe1 = Fe1 + gw(jj) * R' * P * nx * R * det_Jacobi;
%         Fe2 = Fe2 + gw(jj) * R' * P * ny * R * det_Jacobi;
        for aa = 1: 4
            Na = R(aa);
            Fe1(aa, 1) = Fe1(aa, 1) + gw(jj) * Na * (R' * P) * nx * det_Jacobi;
            Fe2(aa, 1) = Fe2(aa, 1) + gw(jj) * Na * (R' * P) * ny * det_Jacobi;
        end
    end
end
%%%%%%% 当压力施加于单元第4边时，Fe1和Fe2的计算
%%%%%%% Fe1, Fe2的数值积分