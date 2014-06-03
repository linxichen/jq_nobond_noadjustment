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
zbar xxibar

yss nss kss mmuss css wss;


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

model(linear);
	//1. HH Budget
	nss*w + wss*n + d = c;

	//2. Labor Demand
	w = -(1-ttheta)*yss/nss*mmu + (1-mmuss)*(1-ttheta)/n*y - nss^(-2)*(1-mmuss)*(1-ttheta)*yss*n;

	//3. Labor Supply
	w = aalpha/(1-nss)*c + aalpha*c/((1-nss)^2)*n;

	//4. Capital Demand
	-css^(-2)*(1-mmuss*xxibar)*c - xxibar/css*mmu - mmuss*xxibar/css*xxi = bbeta*mk(+1);

	// 5. Resource Constraint 
	(1-ddelta)*k(-1) + y - k - c = 0;

	// 6. Financial Constraint
	xxibar*kss*xxi + xxibar*k =  y;

	// 7. Output Definiton
	y = yss*z + ttheta*yss/kss*k(-1) + (1-ttheta)*yss/nss*n;

	// 8. Investment
	inv = k-(1-ddelta)*k(-1);

	// 9. MK Definition
	mk = -(css^(-2))*(1-ddelta+(1-mmuss)*ttheta*yss/kss)*c - ttheta*yss/(css*kss)*mmu + (1-mmuss)*ttheta/(css*kss)*y - (1-mmuss)*ttheta*yss/(css*kss^2)*k(-1);

	// 10 LOM 1
	z = rrhozz*z(-1) + rrhozxxi*xxi(-1) + ssigmaepsz*epsz;

	// 11 LOM 2
	xxi = rrhoxxiz*z(-1) + rrhoxxixxi*xxi(-1) + ssigmaepsxxi*epsxxi;

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
save
// stoch_simul(hp_filter=1600,order = 1,periods=100000,irf=40,drop=2000,loglinear) y c inv n w z d; % compute polices up to 1st order