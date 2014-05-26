function f = case2(n)
    k = 6;
    m = 3;
    ddelta = 0.025;
    aalpha = -.7;
    ttheta = 0.36;
    xxi = 0.07;
    z = 1.1;
    
    f = (1-ddelta)*k*m - (1-ddelta)*(k*m+aalpha*ttheta/(1-ttheta))*n ...
      + (1-1/xxi)*z*k^ttheta*m*n^(1-ttheta) ...
      + (1/xxi-1)*k^ttheta*(m + aalpha*ttheta/((1-ttheta)*k))*n^(2-ttheta) ...
      - 1 + ddelta + (1-ddelta)*n;
end