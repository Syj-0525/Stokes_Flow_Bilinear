function [Be1, Be2] = LocAssem_Be(eleCrtlPts)
%%%%%%% 初始化Be1和Be2
Be1 = zeros(4, 4);
Be2 = zeros(4, 4);
%%%%%%% 初始化Be1和Be2

%%%%%%% 高斯积分点和权重
nqp = 10;
[gp, gw] = Gauss(nqp, -1, 1);
% gp = [0.9324695142031521, 0.6612093864662645, 0.2386191860831969, -0.9324695142031521, -0.6612093864662645, -0.2386191860831969];
% gw = [0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.1713244923791704, 0.3607615730481386, 0.4679139345726910];
kesi = gp;
ita = gp;
%%%%%%% 高斯积分点和权重 

%%%%%%% Be1和Be2的数值积分
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
                    1/4 * (1 + kesi(ii))
                    1/4 * (1 - kesi(ii)) ];
        %%%%%%% 速度插值函数对kesi和ita导数

        %%%%%% Jacobi相关计算
        dx_dkesi = R_kesi' * eleCrtlPts(:,1);
        dx_dita = R_ita' * eleCrtlPts(:,1);
        dy_dkesi = R_kesi' * eleCrtlPts(:,2);
        dy_dita = R_ita' * eleCrtlPts(:,2);
        Jacobi = [ dx_dkesi dy_dkesi
                   dx_dita  dy_dita ];
        R_xy = inv(Jacobi) * [R_kesi'; R_ita'];
        R_x = R_xy(1,:)';
        R_y = R_xy(2,:)';
        det_Jacobi = det(Jacobi);
        %%%%%% Jacobi相关计算

        %%%%%% Be1和Be2单元子块矩阵计算
%         Be1 = Be1 + gw(ii) * gw(jj)  * R_x * R' * det_Jacobi;
%         Be2 = Be2 + gw(ii) * gw(jj)  * R_y * R' * det_Jacobi;
        for aa = 1: 4
            Na_x = R_x(aa); Na_y = R_y(aa);
            for bb = 1: 4
                Nb = R(bb);
                Be1(aa, bb) = Be1(aa, bb) + gw(ii) * gw(jj) * Na_x * Nb * det_Jacobi;
                Be2(aa, bb) = Be2(aa, bb) + gw(ii) * gw(jj) * Na_y * Nb * det_Jacobi;
            end
        end
        %%%%%% Be1和Be2单元子块矩阵计算
    end
end
%%%%%%% Be1和Be2的数值积分