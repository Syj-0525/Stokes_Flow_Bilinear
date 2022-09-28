function [Ce1, Ce2] = LocAssem_Ce(eleCtrlPts)
%%%%%%% 初始化Ce1和Ce2
Ce1 = zeros(4, 4);
Ce2 = zeros(4, 4);
%%%%%%% 初始化Ce1和Ce2

%%%%%%% 高斯积分点和权重
nqp = 10;
[gp, gw] = Gauss(nqp, -1, 1);
% gp = [0.9324695142031521, 0.6612093864662645, 0.2386191860831969, -0.9324695142031521, -0.6612093864662645, -0.2386191860831969];
% gw = [0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.1713244923791704, 0.3607615730481386, 0.4679139345726910];
kesi = gp;
ita = gp;
%%%%%%% 高斯积分点和权重

%%%%%%% Ce1和Ce2的数值积分
for ii = 1: nqp
    for jj = 1: nqp
        %%%%%%% 形函数
        R = [1/4 * (1 - kesi(ii)) * (1 - ita(jj))
             1/4 * (1 + kesi(ii)) * (1 - ita(jj))
             1/4 * (1 + kesi(ii)) * (1 + ita(jj))
             1/4 * (1 - kesi(ii)) * (1 + ita(jj))];
        %%%%%%% 形函数

        %%%%%%% 形函数对kesi和ita导数
        R_kesi = [ - 1/4 * (1 - ita(jj))
                     1/4 * (1 - ita(jj)) 
                     1/4 * (1 + ita(jj))
                   - 1/4 * (1 + ita(jj)) ];
        R_ita = [ - 1/4 * (1 - kesi(ii))
                  - 1/4 * (1 + kesi(ii))
                    1/4 * (1 + kesi(jj))
                    1/4 * (1 - kesi(jj)) ];
        %%%%%%% 速度插值函数对kesi和ita导数

        %%%%%% Jacobi相关计算
        dx_dkesi = R_kesi' * eleCtrlPts(:,1);
        dx_dita = R_ita' * eleCtrlPts(:,1);
        dy_dkesi = R_kesi' * eleCtrlPts(:,2);
        dy_dita = R_ita' * eleCtrlPts(:,2);
        Jacobi = [ dx_dkesi dy_dkesi
                   dx_dita  dy_dita];
        R_xy = inv(Jacobi) * [R_kesi'; R_ita'];
        R_x = R_xy(1,:)';
        R_y = R_xy(2,:)';
        det_Jacobi = det(Jacobi);
        %%%%%% Jacobi相关计算

        %%%%%% Ce1和Ce2单元子块矩阵计算
%         Ce1 = Ce1 + gw(ii) * gw(jj)  * R * R_x' * det_Jacobi;
%         Ce2 = Ce2 + gw(ii) * gw(jj)  * R * R_y' * det_Jacobi;
        for aa = 1: 4
            Na = R(aa);
            for bb = 1: 4
                Nb_x = R_x(bb); Nb_y = R_y(bb);
                Ce1(aa, bb) = Ce1(aa, bb) + gw(ii) * gw(jj) * Na * Nb_x * det_Jacobi;
                Ce2(aa, bb) = Ce2(aa, bb) + gw(ii) * gw(jj) * Na * Nb_y * det_Jacobi;
            end
        end
        %%%%%% Ce1和Ce2单元子块矩阵计算
    end
end
%%%%%%% Ce1和Ce2的数值积分