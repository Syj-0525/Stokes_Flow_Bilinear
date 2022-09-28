function [Ae11, Ae12, Ae21, Ae22] = LocAssem_Ae(eleCtrlPts, viscosity)
%%%%%%% 初始化Ae11, Ae12, Ae21和Ae22
Ae11 = zeros(4, 4);
Ae12 = zeros(4, 4);
Ae21 = zeros(4, 4);
Ae22 = zeros(4, 4);
%%%%%%% 初始化Ae11, Ae12, Ae21和Ae22

%%%%%%% 高斯积分点和权重
nqp = 10;
[gp, gw] = Gauss(nqp, -1, 1);
% gp = [0.9324695142031521, 0.6612093864662645, 0.2386191860831969, -0.9324695142031521, -0.6612093864662645, -0.2386191860831969];
% gw = [0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.1713244923791704, 0.3607615730481386, 0.4679139345726910];
kesi = gp;
ita = gp;
%%%%%%% 高斯积分点和权重

%%%%%%% Ae11, Ae12, Ae21和Ae22的数值积分
for ii = 1: nqp
    for jj = 1: nqp
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

        %%%%%% Ae11, Ae12, Ae21和Ae22单元子块矩阵计算
%         Ae11 = Ae11 + viscosity * gw(ii) * gw(jj) * ( (R_x * R_x') + R_y * R_y') * det_Jacobi;
%         Ae22 = Ae22 + viscosity * gw(ii) * gw(jj) * ( (R_y * R_y') + R_x * R_x') * det_Jacobi;  
        for aa = 1: 4
            Na_x = R_x(aa); Na_y = R_y(aa);
            for bb = 1: 4
                Nb_x = R_x(bb); Nb_y = R_y(bb);
                Ae11(aa, bb) = Ae11(aa, bb) + viscosity * gw(ii) * gw(jj) * ( (Na_x * Nb_x) + (Na_y * Nb_y) ) * det_Jacobi; 
                Ae22(aa, bb) = Ae22(aa, bb) + viscosity * gw(ii) * gw(jj) * ( (Na_x * Nb_x) + (Na_y * Nb_y) ) * det_Jacobi; 
            end
        end
        %%%%%% Ae11, Ae12, Ae21和Ae22单元子块矩阵计算
    end
end
%%%%%%% Ae11, Ae12, Ae21和Ae22的数值积分