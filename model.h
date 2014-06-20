#include <iostream>
#include <iomanip>
#include <fstream>

// Define an class that contains parameters and steady states
struct para {
	// Model parameters
	double aalpha;
	double bbeta ;
	double ddelta;
	double ttheta;
	double kkappa;
	double ttau  ;
	double xxibar;
	double zbar  ;
	double dbar  ;
	double rrhozz;
	double rrhozxxi;
	double rrhoxxiz;
	double rrhoxxixxi;
	double var_epsz;
	double var_epsxxi;
	double A[4];
	double Ssigma_e[4];

	// Steady States
	double kss;
	double nss;
	double css;
	double wss;
	double dss;
	double mmuss;
	double mkss;
	double yss;


	// Find steady state and find aalpha based steady state target
	__host__ __device__
	void complete() {
		// Fill A and Ssigma_e
		A[0] = rrhozz; A[2] = rrhozxxi;
		A[1] = rrhoxxiz; A[3] = rrhoxxixxi;
		Ssigma_e[0] = var_epsz*var_epsz;
		Ssigma_e[3] = var_epsxxi*var_epsxxi;

		// Find aalpha based on SS computation
		double kovern = pow(xxibar,1/(ttheta-1));
		double covern = pow(kovern,ttheta) - ddelta*kovern;
		mmuss = 1 - ( bbeta*(1-ddelta)-1+xxibar )/( xxibar*(1-bbeta*ttheta)  );
		aalpha = double(0.7/0.3)*(1/covern)*(1-mmuss)*(1-ttheta)*pow(kovern,ttheta);
		double G = ( (1-mmuss)*(1-ttheta)*pow(kovern,ttheta) ) / ( aalpha*covern );
		nss = G/(1+G);
		css = nss*covern;
		kss = nss*kovern;
		wss = aalpha*css/(1-nss);
		dss = css - wss*nss;
		mkss = (1-ddelta+(1-mmuss)*zbar*ttheta*pow(kss,ttheta-1)*pow(nss,1-ttheta))/css;
		yss = zbar*pow(kss,ttheta)*pow(nss,1-ttheta);
	};

	// Export parameters to a .m file in MATLAB syntax
	__host__
	void exportmatlab(std::string filename) {
		std::ofstream fileout(filename.c_str(), std::ofstream::trunc);

		// Model Parameters
		fileout << std::setprecision(16) << "aalpha=" << aalpha << ";"<< std::endl;
		fileout << std::setprecision(16) << "bbeta=" << bbeta << ";"<< std::endl;
		fileout << std::setprecision(16) << "ddelta=" << ddelta << ";"<< std::endl;
		fileout << std::setprecision(16) << "ttheta=" << ttheta << ";"<< std::endl;
		fileout << std::setprecision(16) << "xxibar=" << xxibar << ";"<< std::endl;
		fileout << std::setprecision(16) << "zbar=" << zbar << ";"<< std::endl;
		fileout << std::setprecision(16) << "rrhozz=" << rrhozz << ";"<< std::endl;
		fileout << std::setprecision(16) << "rrhozxxi=" << rrhozxxi << ";"<< std::endl;
		fileout << std::setprecision(16) << "rrhoxxiz=" << rrhoxxiz << ";"<< std::endl;
		fileout << std::setprecision(16) << "rrhoxxixxi=" << rrhoxxixxi << ";"<< std::endl;
		fileout << std::setprecision(16) << "ssigmaepsz=" << sqrt(var_epsz) << ";"<< std::endl;
		fileout << std::setprecision(16) << "ssigmaepsxxi=" << sqrt(var_epsxxi) << ";"<< std::endl;

		// Steady States
		fileout << std::setprecision(16) << "kss=" << kss << ";"<< std::endl;
		fileout << std::setprecision(16) << "nss=" << nss << ";"<< std::endl;
		fileout << std::setprecision(16) << "css=" << css << ";"<< std::endl;
		fileout << std::setprecision(16) << "wss=" << wss << ";"<< std::endl;
		fileout << std::setprecision(16) << "dss=" << dss << ";"<< std::endl;
		fileout << std::setprecision(16) << "mmuss=" << mmuss << ";"<< std::endl;
		fileout << std::setprecision(16) << "mkss=" << mkss << ";"<< std::endl;
		fileout << std::setprecision(16) << "yss=" << yss << ";"<< std::endl;
		fileout.close();
	};
};

void linearizedmodel(double* A, double* B, double* C, double* rrho, int n, int n_shock, para p) {
	// HH Budget. Correct.
	B[0+3*n] = p.nss;
	B[0+2*n] = p.wss;
	B[0+4*n] = 1;
	B[0+1*n] = -1;

	// Labor Demand. Correct
	B[1+5*n] = (p.ttheta-1)*p.yss/p.nss;
	B[1+6*n] = (1-p.ttheta)*(1-p.mmuss)/p.nss;
	B[1+2*n] = -(1-p.ttheta)*(1-p.mmuss)*p.yss/(p.nss*p.nss);
	B[1+3*n] = -1;

	// Labor Supply. Correct
	B[2+1*n] = p.aalpha/(1-p.nss);
	B[2+2*n] = p.aalpha*p.css/((1-p.nss)*(1-p.nss));
	B[2+3*n] = -1;

	// Capital Demand. Correct.
	A[3+8*n] = p.bbeta; 
	B[3+1*n] = -(1-p.mmuss*p.xxibar)/(p.css*p.css); 
	B[3+5*n] = -p.xxibar/p.css; 
	C[3+1*n] = -p.mmuss*p.xxibar/p.css;

	// Resource Constraint. Correct
	A[4+0*n] = 1; 
	B[4+0*n] = 1-p.ddelta; 
	B[4+6*n] = 1; 
	B[4+1*n] = -1;

	// Financial Constraint. Fixed.
	A[5+0*n] = p.xxibar;
	B[5+6*n] = 1;
	C[5+1*n] = -p.xxibar*p.kss;

	// Output Definition. Correct
	C[6+0*n] = p.yss;
	B[6+0*n] = p.ttheta*p.yss/p.kss;
	B[6+2*n] = (1-p.ttheta)*p.yss/p.nss;
	B[6+6*n] = -1;

	// Investment Definition. Correct
	A[7+0*n] = 1;
	B[7+7*n] = 1;
	B[7+0*n] = 1-p.ddelta;

	// MK defintion:
	B[8+1*n] = -pow(p.css,-2)*(1-p.ddelta+(1-p.mmuss)*p.ttheta*p.yss/p.kss); 
	B[8+5*n] = -p.ttheta*p.yss/(p.css*p.kss); 
	B[8+6*n] = (1-p.mmuss)*p.ttheta/(p.css*p.kss); 
	B[8+0*n] = -(1-p.mmuss)*p.ttheta*p.yss*pow(p.kss,-2)/p.css;
	B[8+8*n] = -1;

	for (int i=0; i< n_shock*n_shock; i++) {
		rrho[i] = p.A[i];
	};
};

// Define state struct that contains "natural" state 
struct state {
	// Data member
	double k, z, xxi, zkttheta;

	// Constructor
	__host__ __device__
	state(double _k, double _z, double _xxi, para p) {
		k = _k;
		z = _z;
		xxi = _xxi;
		zkttheta = _z*pow(_k,p.ttheta);
	};

	// Alternative constructor
	__host__ __device__
	state(double _k, double _z, double _xxi, double _zkttheta) {
		k = _k;
		z = _z;
		xxi = _xxi;
		zkttheta = _zkttheta;
	};
	__host__ __device__
	void load(double _k, double _z, double _xxi, para p) {
		k = _k;
		z = _z;
		xxi = _xxi;
		zkttheta = _z*pow(_k,p.ttheta);
	};
};

// Define shadow struct that contains all shadow values
struct shadow {
	// Data member
	double m1;

	// Constructor
	__host__ __device__
	shadow(double _m1) {
		m1 = _m1;
	};

	__host__ __device__
	void load(double _m1) {
		m1 = _m1;
	};
};

struct binding_hour {
	// Data Member are const coefficents and some model parameters
	double c0, c1, c_oneminusttheta, c_twominusttheta;
	double ttheta;

	// Constructot (Adrian/Projection)
	__host__ __device__
	binding_hour(state s, shadow m, para _p) {
		ttheta = _p.ttheta;
		double aalpha = _p.aalpha;
		double ddelta = _p.ddelta;
		c0 = (1-ddelta)*s.k*m.m1 - 1 + ddelta;
		c1 = (1-ddelta)*(1-s.k*m.m1-aalpha*ttheta/(1-ttheta));
		c_oneminusttheta = (1-1/s.xxi)*m.m1*s.zkttheta;
		c_twominusttheta = (1/s.xxi-1)*s.zkttheta*(m.m1+aalpha*ttheta/(s.k*(1-ttheta)));
	};

	// The function of hour
	__host__ __device__
	double operator()(double n) {
		return c0 + c1*n + c_oneminusttheta*pow(n,1-ttheta) + c_twominusttheta*pow(n,2-ttheta);
	};

	__host__ __device__
	// The derivative of function
	double prime(double n) {
		return c1 + (1-ttheta)*c_oneminusttheta*pow(n,-ttheta) + (2-ttheta)*c_twominusttheta*pow(n,1-ttheta);
	};
};

struct notbinding_hour {
	// Data Member are const coefficents and some model parameters
	double c0,  c_oneminusttheta, c_minusttheta;
	double ttheta;

	// Constructor (Adrian's method/Projection)
	__host__ __device__
	notbinding_hour(state s, shadow m, para _p) {
		double ddelta = _p.ddelta;
		ttheta = _p.ttheta;
		double aalpha = _p.aalpha;
		c0 = (1-ddelta)/m.m1;
		c_oneminusttheta = s.zkttheta*(ttheta/(m.m1*s.k)+(1-ttheta)/aalpha);
		c_minusttheta = -(1-ttheta)*s.zkttheta/aalpha;
	};

	// Constructor (Value Function Iteration)
	__host__ __device__
	notbinding_hour(state s, state splus, para _p) {
		ttheta = _p.ttheta;
		double ddelta = _p.ddelta;
		double aalpha = _p.aalpha;
		c0 = aalpha*(1-ddelta)*s.k - aalpha*splus.k;
		c_oneminusttheta = (aalpha+1-ttheta)*s.zkttheta;
		c_minusttheta = (ttheta-1)*s.zkttheta;
	};

	// The function of hour
	__host__ __device__
	double operator()(double n) {
		return c0 + c_oneminusttheta*pow(n,1-ttheta) + c_minusttheta*pow(n,-ttheta);
	};

	__host__ __device__
	// The derivative of function
	double prime(double n) {
		return (1-ttheta)*c_oneminusttheta*pow(n,-ttheta) + (-ttheta)*c_minusttheta*pow(n,-ttheta-1);
	};
};

// struct of variables implied by natural state, optionally shadow value/kplus tomorrow
struct control {
	// Data member
	double kplus, c, n, w, d, mmu, Y, lhs1;

	// finding the control variables (adrian's method)
	__host__ __device__
	void compute(state s, shadow m, para p, int binding) {
		if (binding == 1) {
			// Case 1: Binding
			n = newton(binding_hour(s,m,p),1e-5,1.0-1e-5,0.3);
			Y = s.zkttheta*pow(n,1-p.ttheta);
			double MPK = p.ttheta*Y/s.k;
			kplus = Y/s.xxi;
			c = Y+(1-p.ddelta)*s.k-kplus;
			mmu = 1-(m.m1*c-1+p.ddelta)/MPK;	
			w = (1-mmu)*(1-p.ttheta)*Y/n;
			lhs1 = (1-mmu*s.xxi)/c;
			d = c - w*n;
		};

		if (binding == 0) {
			// Case 2: Not Binding
			n = newton(notbinding_hour(s,m,p),1e-5,1.0-1e-5,0.3);
			Y = s.zkttheta*pow(n,1-p.ttheta);
			double MPK = p.ttheta*Y/s.k;
			mmu = 0;
			c = (1-p.ddelta+MPK)/m.m1;	
			kplus = (1-p.ddelta)*s.k + Y - c;
			w = (1-p.ttheta)*Y/n;
			lhs1 = (1-mmu*s.xxi)/c;
			d = c - w*n;
		};
	};
	
	// finding the control variables (VFI)
	__host__ __device__
	void compute(state s, state splus, para p, int binding) {
		if (binding == 1) {
			// Case 1: Binding
			double k = s.k;
			double xxi = s.xxi;
			double zkttheta = s.zkttheta;
			kplus = splus.k;
			n = pow(xxi*kplus/zkttheta,1/(1-p.ttheta));
			Y = xxi*kplus;
			c = Y + (1-p.ddelta)*k - kplus;
			mmu = 1-p.aalpha*c/( (1-n)*(1-p.ttheta)*Y/n );
		};

		if (binding == 0) {
			// Case 2: Not Binding
			double k = s.k;
			double zkttheta = s.zkttheta;
			kplus = splus.k;
			n = newton(notbinding_hour(s,splus,p),1e-5,1.0-1e-5,0.3);
			Y = zkttheta*pow(n,1-p.ttheta);
			mmu = 0;
			c = Y + (1-p.ddelta)*k - kplus; 
		};
	};
};
