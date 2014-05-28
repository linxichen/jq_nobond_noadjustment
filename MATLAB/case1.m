function f = case1(n)
    k = 10.955;
    m = 1.2259;
    ddelta = 0.025;
    aalpha = 1.7849;
    ttheta = 0.36;
    xxi = 0.1;
    z = 1;
    
    f = (1-ddelta)*k*m -1 + ddelta...
      + (1-ddelta)*(1 - k*m -aalpha*ttheta/(1-ttheta))*n ...
      + (1-1/xxi)*z*k^ttheta*m*n^(1-ttheta) ...
      + (1/xxi-1)*z*k^ttheta*(m + aalpha*ttheta/((1-ttheta)*k))*n^(2-ttheta);
end