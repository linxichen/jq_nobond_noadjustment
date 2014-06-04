//----------------------------------------------------------------
// 1. Declare variables
//----------------------------------------------------------------
var 

// endogenous variables
w           // variable 1
c           // variable 2           
n           // variable 3
d           // variable 4
k           // variable 5
mmu         // variable 6
y           // variable 7
inv         // variable 8
mk			// variable 9
yovern

// exogenous variables
z			// variable 14
xxi;        // variable 15




//----------------------------------------------------------------
// 2. Exogenous shocks
//----------------------------------------------------------------

varexo epsz epsxxi;

//----------------------------------------------------------------
// 3. Parameters
//----------------------------------------------------------------

parameters 

bbeta aalpha ttheta ddelta
rrhozz rrhoxxixxi rrhozxxi rrhoxxiz
ssigmaepsz ssigmaepsxxi 
zbar xxibar;


//----------------------------------------------------------------
// 4. Calibration 
//----------------------------------------------------------------
cd ../MATLAB;
mypara;
cd ../Dynare
for i=1:length(M_.params)
    deep_parameter_name = M_.param_names(i,:);
    eval(['M_.params(i)  = ' deep_parameter_name ' ;'])
end
 
//----------------------------------------------------------------
// 5. Steady State
//----------------------------------------------------------------
steady_state_model;
    xxi = xxibar;
    z = zbar;
    kovern = xxi^(1/(ttheta-1));
    covern = kovern^ttheta - ddelta*kovern;
    mmu = 1 - (bbeta*(1-ddelta)-1+xxi)/(xxi*(1-bbeta*ttheta));
    G = (1-mmu)*(1-ttheta)*kovern^ttheta / (aalpha*covern);
    n = G/(1+G);
    c = covern*n;
	w = aalpha*c/(1-n);
    k = kovern*n;
    d = c - w*n;
	inv = ddelta*k;
    mk = (c^-1)*( (1-ddelta) + (1-mmu)*ttheta*z*k^(ttheta-1)*n^(1-ttheta) ) ;
	y = z*k^ttheta*n^(1-ttheta);
    yovern = y/n;
end;
steady;


//----------------------------------------------------------------
// 6. Model
//----------------------------------------------------------------

model;
	//1. HH Budget
	w*n + d - c = 0;

	//2. Labor Demand
	(1-mmu)*(1-ttheta)*z*(k(-1)^ttheta)*n^(-ttheta) = w;

	//3. Labor Supply
	w = aalpha*c^(1)/(1-n);

	//4. Capital Demand
	( 1-mmu*xxi )*c^(-1)
	= bbeta*(c(+1)^(-1))*( 1-ddelta+(1-mmu(+1))*z(+1)*ttheta*k^(ttheta-1)*n(+1)^(1-ttheta) );

	// 5. Resource Constraint 
	(1-ddelta)*k(-1) + y - k - c = 0;

	// 6. Financial Constraint
	xxi*( k ) =  y;

	// 7. Output
	y = z*(k(-1))^(ttheta)*n^(1-ttheta);

	// 8. Investment
	inv = k-(1-ddelta)*k(-1);

	// 9. MK Definition
	mk = (c^(-1))*( (1-ddelta)+(1-mmu)*z*ttheta*k(-1)^(ttheta-1)*n^(1-ttheta) );

	// 10
	(z/zbar) = (z(-1)/zbar)^(rrhozz)*(xxi(-1)/xxibar)^(rrhozxxi)*exp(ssigmaepsz*epsz);

	// 11
	(xxi/xxibar)= (xxi(-1)/xxibar)^(rrhoxxixxi)*(z(-1)/zbar)^(rrhoxxiz)*exp(ssigmaepsxxi*epsxxi);

    // 12
    yovern = y/n;
end;

//----------------------------------------------------------------
// 7. Computation
//----------------------------------------------------------------

/*****
initval;
w = (wss);  
c = (css);       
// n = (nss);          
// R = (Rss);      
// b = (bss);
d = (dss);  
k = (kss);  
inv = (invss);
y = (yss);  
mmu = (mmuss);
z = (zss);
xxi = (xxiss);
mk = (mkss);
// mb = (mbss);
// mc = (mcss);
epsz = 0;
epsxxi = 0;
end;
*****/

shocks;
var epsz;
stderr 1;
var epsxxi;
stderr 1;
end;

stoch_simul(order = 1,periods=100000,irf=40,drop=2000, hp_filter=1600) y c inv n w mmu z d; % compute polices up to 1st order
// stoch_simul(hp_filter=1600,order = 1,periods=100000,irf=40,drop=2000,loglinear) y c inv n w z d; % compute polices up to 1st order