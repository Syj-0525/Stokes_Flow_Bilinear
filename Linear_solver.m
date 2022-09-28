function [x] = Linear_solver(K, B)
x = K\B;   % 高斯消元