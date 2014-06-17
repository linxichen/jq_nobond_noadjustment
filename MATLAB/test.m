clear all
close all
clc

%%
coeff = ones(3,3);
coeff(1) = 1;
coeff(2) = 2;
coeff(3) = 3;
coeff(4) = 4;
coeff(5) = 5;
coeff(6) = 6;
coeff(7) = 7;
coeff(8) = 8;
coeff(9) = 9;

eval = 0;
for i_k = 1:3
    for i_z = 1:3
        eval = eval + coeff(i_k,i_z)*polyval(ChebyshevPoly(i_k-1),0.3)*polyval(ChebyshevPoly(i_z-1),0.4);
    end
end