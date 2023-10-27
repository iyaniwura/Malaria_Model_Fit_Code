// applicable to weekly cases data by symptoms onset date
functions {
	real[] SEIR(real t, 
							real[] y, 
							real[] theta, 
							real[] x_r,
							int[] x_i
	) {
		int TotalPop = x_i[1];            // total population
		
		real L_h = x_r[1];               //  
		real beta_h = x_r[2];               //  
		real mu_h = x_r[3];          // 
		real sigma_h = x_r[4];          // 
		real gamma_h = x_r[5];          // 
		real beta_v = x_r[6];          // 
		real mu_v = x_r[7];          // 
		real PIE = x_r[8];

		real y0_vars[7];

		row_vector[7] dydt;       //  an array containing the derivative of all the compartment arrange in the order S, E1, E2, I1, I2, R, and Inci 
		
		int Nh;
		real b;
		real lambda_h;
		real lambda_v;
		
		real Sh = y[1]; 
		real Eh = y[2];
		real Ih = y[3]; 
		real Rh = y[4]; 
		real Sv =  y[5]; 
		real Iv =  y[6];
		real CumInci = y[7];
		
		real dSh_dt;
		real dEh_dt;
		real dIh_dt;
		real dRh_dt;
		real dSv_dt;
		real dIv_dt;
		real dInci_dt;
		
		real i0 = theta[1];
		real epsilon = theta[2];
		real p = theta[3];
		real b0 = theta[4];
		real Iv0 = theta[5];
		real w = theta[6];
		real w_h = theta[7];
		real delta = theta[8];
		real L_v = theta[9];
		real tc = theta[10];

		y0_vars = x_r[9:15];
		
		
		Nh = TotalPop;
		b = b0*(1 + epsilon*cos(2*PIE*w*(t + tc) ) );
		lambda_h = b*beta_h*Iv/Nh;
		lambda_v = b*beta_v*(Ih)/Nh;
		
		dSh_dt = L_h - lambda_h*Sh + w_h*Rh - mu_h*Sh;    
		dEh_dt = lambda_h*Sh - sigma_h*Eh - mu_h*Eh;
		dIh_dt = sigma_h*Eh - gamma_h*Ih - (mu_h + delta)*Ih;
		dRh_dt = gamma_h*Ih - w_h*Rh - mu_h*Rh;
		dSv_dt = L_v - lambda_v*Sv - mu_v*Sv;
		dIv_dt = lambda_v*Sv - mu_v*Iv;
		dInci_dt = sigma_h*Eh;
		
		return {dSh_dt, dEh_dt, dIh_dt, dRh_dt, dSv_dt, dIv_dt, dInci_dt};
		
	}
}



data {
	// Structure
	int TotalPop;
	int n_months;       // total number of weeks
	real t0;            //starting time
	real ts[n_months];  
	int doprint;
	int inference;
	real L_h;
	real beta_h;
	real mu_h;
	real sigma_h;
	real gamma_h;
	real beta_v;
	real mu_v;
	real PIE;

	
	
	// Data to fit
	int Reported_cases[n_months];
	
	// Priors
	real p_i0[2];
	real p_epsilon[2];
	real p_pp[2];
	real p_b0[2];
	real p_Iv0[2];
	real p_w[2];
	real p_wh[2];
	real p_delta[2];
	real p_Lv[2];
	real p_tc[2];
	real p_phi;

	// Fixed parameters
	  real y0_vars[7];  
}
	

transformed data {
	real x_r[15];
	int x_i[1];
	
	x_i[1] = TotalPop; 
	
	x_r[1] =  L_h;
	x_r[2] =  beta_h;
	x_r[3] =  mu_h;
	x_r[4] =  sigma_h;
	x_r[5] =  gamma_h;
	x_r[6] =  beta_v;
	x_r[7] =  mu_v;
	x_r[8] =  PIE;
	
	x_r[9:15] = y0_vars;
	
}



parameters{
	real<lower=1, upper=5000> i0;  
	real<lower=0.001, upper=5> epsilon;  
	real<lower=0, upper=1> p;  
	real<lower=1, upper=500> b0; 
	real<lower=1, upper=50000> Iv0; 
	real<lower=0, upper=2> w; 
	real<lower=0, upper=2> w_h; 
	real<lower=0, upper=2> delta; 
	real<lower=0, upper=50000> L_v; 
	real<lower=-6.283185, upper=2*6.283185> tc; 
	real<lower=0> phi; // 
}

transformed parameters {
	// transformed parameters
	real y[n_months, 7]; // raw ODE output
	real y0_temp[7];     // initial state
	real y0[7];  //
	
	// change of format for integrate_ode_rk45
	real theta[10];     // vector of parameters
	
	
	// ode outputs
	// cumulative incidence
	vector[n_months] Cum_inci; // overall case incidence by month
	
	// incidence
	vector[n_months] Inci; // overall case incidence by month
	
		
	// set up the initial conditions:
	y0_temp[1] = y0_vars[1] - (y0_vars[2]*i0 + y0_vars[3]*i0);
	y0_temp[2] = y0_vars[2]*i0;
	y0_temp[3] = y0_vars[3]*i0;
	y0_temp[4] = y0_vars[4];
	y0_temp[5] = y0_vars[5];
	y0_temp[6] = Iv0;  
	y0_temp[7] = y0_vars[7];

	
	y0 = to_array_1d(y0_temp);
	
	// change of format for integrate_ode_rk45
	theta[1:10] =  { i0, epsilon, p, b0, Iv0, w, w_h, delta, L_v,  tc };

	// // // run ODE solver
	y = integrate_ode_rk45( SEIR, // ODE function
													y0, // initial states
													t0, // t0
													ts, // evaluation dates (ts)
													theta, // parameters
													x_r, // real data
													x_i, // integer data
													1.0E-6, 1.0E-6, 1.0E3); // tolerances and maximum steps
	// y = integrate_ode_bdf(SEIR, y0, t0, ts, theta, x_r, x_i); // tolerances and maximum steps

	// Extracting incidence from the output of the ode solver	
	Cum_inci = to_vector(y[,7]) ;  
	
	Inci[1] =  Cum_inci[1]; 	  Inci[2:n_months] =  Cum_inci[2:n_months] - Cum_inci[1:(n_months-1)];
	
}

model {
	
	// priors
	i0 ~ lognormal(p_i0[1], p_i0[2]);
	epsilon ~ lognormal(p_epsilon[1], p_epsilon[2]);
	p ~ normal(p_pp[1], p_pp[2]);
	b0 ~ lognormal(p_b0[1], p_b0[2]);
	Iv0 ~ lognormal(p_Iv0[1], p_Iv0[2]);
	w ~ lognormal(p_w[1], p_w[2]);
	w_h ~ lognormal(p_wh[1], p_wh[2]);
	delta ~ lognormal(p_delta[1], p_delta[2]);
	L_v ~ lognormal(p_Lv[1], p_Lv[2]);
	tc ~ normal(p_tc[1], p_tc[2]);
	phi ~ exponential(p_phi);


	if (inference!=0){
			Reported_cases ~ neg_binomial_2(Inci*(1 - p), 1/phi);
			
		}
}


generated quantities{
	real pred_cases[n_months];
	pred_cases = neg_binomial_2_rng(Inci*(1 - p), 1/phi);

}
