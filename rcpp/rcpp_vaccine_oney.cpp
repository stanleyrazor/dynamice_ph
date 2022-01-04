#define _USE_MATH_DEFINES
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_vaccine_oney(NumericMatrix in_Comp, List parm, List siaparm, NumericMatrix beta_full,
	NumericVector pop_full, NumericVector cov1, double cov2, int t_start)
{
	/* This function runs measles transmission, ageing, and vaccination at a specific age group and timestep
	during calendar years for sumulation. */

	// =================================================
	// Set-up (index, parameters)
	// =================================================
	int i_M = 0;    // Maternal Immune
	int i_S = 1;    // Susceptible
	int i_I = 2;    // Infectious
	int i_R = 3;    // Recovered
	int i_V1S = 4;  // Vaccinated susceptible - 1 dose
	int i_V1I = 5;  // Vaccinated infectious  - 1 dose
	int i_V1R = 6;  // Vaccinated recovered   - 1 dose
	int i_V2S = 7;  // Vaccinated susceptible - 2 dose
	int i_V2I = 8;  // Vaccinated infectious  - 2 dose
	int i_V2R = 9;  // Vaccinated recovered   - 2 dose
	int i_V3S = 10; // Vaccinated susceptible - 3 dose
	int i_V3I = 11; // Vaccinated infectious  - 3 dose
	int i_V3R = 12; // Vaccinated recovered   - 3 dose
	int i_V1F = 13; // Vaccinated population  - 1 dose (counter)

	List outp;                           // output list, including compartment distribution and cases
	NumericMatrix trans_Comp(254, 14);   // compartments after including transmission, 254 age groups, 14 states
	NumericMatrix out_Comp(254, 14);     // compartments after including transmission and ageing, 254 age groups, 14 states
	NumericVector betta(254);            // contact rate reported by a contactor of age a, 254 age groups
	NumericVector cyc(254);              // case prevalence (adjusted for seasonality), 254 age groups
	NumericVector newinfect(254);        // new infections/cases, 254 age groups
	double tcycle = 0.0;                 // seasonality
	double lambda = 0.0;                 // force of infection
	double pop_fert_SR = 0.0;            // population of S and R at fertility age
	double pop_fert_R = 0.0;             // population of R at fertility age
	double prp_R = 0.0;                  // proportion of being born with maternal immunity
	double n_eff = 0.0;                  // population with effective protection of MCV1
	double n_v1 = 0.0;                   // population who receive MCV1
	double p_eff = 0.0;                  // proportion of effective vaccine protection among those receive MCV1
	double ve2 = 0.0;                    // 2nd dose vaccine efficay conditioned on MCV1
	
	double tstep = as<double>(parm["tstep"]);                // timesteps per year
	double gamma = as<double>(parm["gamma"]);                // recovery rate per timestep
	double amp = as<double>(parm["amp"]);	                 // amplification for seasonality
	NumericVector ve1 = as<NumericVector>(parm["ve1"]);      // vaccine efficacy of the first dose, 254 age groups
	double ve2plus = as<double>(parm["ve2plus"]);            // vaccine protection for two doses
	double age_w = 52/tstep;                                 // weekly ageing rate per timestep
	double age_y = 1/tstep;  				                 // annualy ageing rate per timestep
	double wane = (12/6)/tstep;  		                     // waning rate maternal immunity per timestep (duration of 6 months)
    int sia_implement = as<int>(siaparm["sia_implement"]);   // SIA implementation methods, 0: no SIA, 1: Portnoy's data, 2: 7.7% never-reach

    for (int t = t_start; t <= (t_start+tstep); ++t)
	{
		// =================================================
        // Transmission
        // =================================================
        tcycle = 1.0 + amp*sin(2.0*M_PI*t/tstep);            //Seasonality in Beta, M_PI = 3.14159265358979323846
        cyc = tcycle*(in_Comp(_,i_I) + in_Comp(_,i_V1I) + in_Comp(_,i_V2I) + in_Comp(_,i_V3I) + 1.0e-9);

		for (int a = 0; a < 254; ++a)
		{
			betta = beta_full(a,_);
            lambda = 1.0 - exp(-sum(betta*cyc));
            newinfect[a] += lambda*(in_Comp(a,i_S) + in_Comp(a,i_V1S) + in_Comp(a,i_V2S) + in_Comp(a,i_V3S));
            //Rcout << "Age group = " << a+1 << "\n" // print out FOI to check
            //      << "FOI = " << lambda << "\n";

            trans_Comp(a,i_M)   = in_Comp(a,i_M);
            trans_Comp(a,i_S)   = in_Comp(a,i_S)   - lambda*in_Comp(a,i_S);
            trans_Comp(a,i_I)   = in_Comp(a,i_I)   + lambda*in_Comp(a,i_S)   - gamma*in_Comp(a,i_I);
            trans_Comp(a,i_R)   = in_Comp(a,i_R)   + gamma *in_Comp(a,i_I);
            trans_Comp(a,i_V1S) = in_Comp(a,i_V1S) - lambda*in_Comp(a,i_V1S);
            trans_Comp(a,i_V1I) = in_Comp(a,i_V1I) + lambda*in_Comp(a,i_V1S) - gamma*in_Comp(a,i_V1I);
            trans_Comp(a,i_V1R) = in_Comp(a,i_V1R) + gamma *in_Comp(a,i_V1I);
            trans_Comp(a,i_V2S) = in_Comp(a,i_V2S) - lambda*in_Comp(a,i_V2S);
            trans_Comp(a,i_V2I) = in_Comp(a,i_V2I) + lambda*in_Comp(a,i_V2S) - gamma*in_Comp(a,i_V2I);
            trans_Comp(a,i_V2R) = in_Comp(a,i_V2R) + gamma *in_Comp(a,i_V2I);
            trans_Comp(a,i_V3S) = in_Comp(a,i_V3S) - lambda*in_Comp(a,i_V3S);
            trans_Comp(a,i_V3I) = in_Comp(a,i_V3I) + lambda*in_Comp(a,i_V3S) - gamma*in_Comp(a,i_V3I);
            trans_Comp(a,i_V3R) = in_Comp(a,i_V3R) + gamma *in_Comp(a,i_V3I);
            trans_Comp(a,i_V1F) = in_Comp(a,i_V1F);
	    }
        //Rcout << "\ntrnansmission cycle completed\n";


        // =================================================
        // Ageing and routine vaccination (MCV1, MCV2)
        // =================================================
        //NumericMatrix out_Comp(254, 14);         // output compartments after including transmission and ageing, 254 age groups, 14 states

        // ageing process: 3-100 years old, no vaccination
        for (int a = 253; a > 155; --a)
		{
			out_Comp(a,i_M)   = trans_Comp(a,i_M)   - age_y*trans_Comp(a,i_M)   + age_y*trans_Comp(a-1,i_M)   - wane*trans_Comp(a,i_M);
            out_Comp(a,i_S)   = trans_Comp(a,i_S)   - age_y*trans_Comp(a,i_S)   + age_y*trans_Comp(a-1,i_S)   + wane*trans_Comp(a,i_M);
            out_Comp(a,i_I)   = trans_Comp(a,i_I)   - age_y*trans_Comp(a,i_I)   + age_y*trans_Comp(a-1,i_I)  ;
            out_Comp(a,i_R)   = trans_Comp(a,i_R)   - age_y*trans_Comp(a,i_R)   + age_y*trans_Comp(a-1,i_R)  ;
            out_Comp(a,i_V1S) = trans_Comp(a,i_V1S) - age_y*trans_Comp(a,i_V1S) + age_y*trans_Comp(a-1,i_V1S);
            out_Comp(a,i_V1I) = trans_Comp(a,i_V1I) - age_y*trans_Comp(a,i_V1I) + age_y*trans_Comp(a-1,i_V1I);
            out_Comp(a,i_V1R) = trans_Comp(a,i_V1R) - age_y*trans_Comp(a,i_V1R) + age_y*trans_Comp(a-1,i_V1R);
            out_Comp(a,i_V2S) = trans_Comp(a,i_V2S) - age_y*trans_Comp(a,i_V2S) + age_y*trans_Comp(a-1,i_V2S);
            out_Comp(a,i_V2I) = trans_Comp(a,i_V2I) - age_y*trans_Comp(a,i_V2I) + age_y*trans_Comp(a-1,i_V2I);
            out_Comp(a,i_V2R) = trans_Comp(a,i_V2R) - age_y*trans_Comp(a,i_V2R) + age_y*trans_Comp(a-1,i_V2R);
            out_Comp(a,i_V3S) = trans_Comp(a,i_V3S) - age_y*trans_Comp(a,i_V3S) + age_y*trans_Comp(a-1,i_V3S);
            out_Comp(a,i_V3I) = trans_Comp(a,i_V3I) - age_y*trans_Comp(a,i_V3I) + age_y*trans_Comp(a-1,i_V3I);
            out_Comp(a,i_V3R) = trans_Comp(a,i_V3R) - age_y*trans_Comp(a,i_V3R) + age_y*trans_Comp(a-1,i_V3R);
            out_Comp(a,i_V1F) = trans_Comp(a,i_V1F) - age_y*trans_Comp(a,i_V1F) + age_y*trans_Comp(a-1,i_V1F);
            //Rcout << a+1 << " ";
        }
        //Rcout << "\nYearly age groups done\n";

        // ageing process and routine vaccination:  2 weeks - 2 years old
        for (int a = 155; a > 0; --a)
		{
			if (a == 71)
			{
				// ageing process and routine vaccination:  72 weeks old
				out_Comp(71,i_M)   = trans_Comp(71,i_M)
								   - age_w*trans_Comp(71,i_M)
								   + age_w*trans_Comp(71-1,i_M)*(1.0-cov1[71-1])
								   - wane *trans_Comp(71,i_M);
				
				out_Comp(71,i_S)   = trans_Comp(71,i_S)
								   - age_w*trans_Comp(71,i_S)
								   + age_w*trans_Comp(71-1,i_S)*(1.0-cov1[71-1])
								   + wane *trans_Comp(71,i_M);
				
				out_Comp(71,i_I)   = trans_Comp(71,i_I)
								   - age_w*trans_Comp(71,i_I)
								   + age_w*trans_Comp(71-1,i_I)*(1.0-cov1[71-1]);
				
				out_Comp(71,i_R)   = trans_Comp(71,i_R)
								   - age_w*trans_Comp(71,i_R)
								    + age_w*trans_Comp(71-1,i_R)*(1.0-cov1[71-1]);
				
				out_Comp(71,i_V1S) = trans_Comp(71,i_V1S)
								- age_w*trans_Comp(71,i_V1S)
								+ age_w*trans_Comp(71-1,i_V1S)*(1.0-cov2)
								+ age_w*(trans_Comp(71-1,i_M)+trans_Comp(71-1,i_S))*cov1[71-1]*(1.0-ve1[71-1]);
				
				out_Comp(71,i_V1I) = trans_Comp(71,i_V1I)
								- age_w*trans_Comp(71,i_V1I)
								+ age_w*trans_Comp(71-1,i_V1I)*(1.0-cov2)
								+ age_w*trans_Comp(71-1,i_I)*cov1[71-1];
				
				out_Comp(71,i_V1R) = trans_Comp(71,i_V1R)
								- age_w*trans_Comp(71,i_V1R)
								+ age_w*trans_Comp(71-1,i_V1R)*(1.0-cov2)
								+ age_w*(trans_Comp(71-1,i_M)+trans_Comp(71-1,i_S))*cov1[71-1]*ve1[71-1]
								+ age_w*trans_Comp(71-1,i_R)*cov1[71-1];
				
				out_Comp(71,i_V1F) = trans_Comp(71,i_V1F)
								- age_w*trans_Comp(71,i_V1F)
								+ age_w*trans_Comp(71-1,i_V1F)*(1.0-cov2)
								+ age_w*(trans_Comp(71-1,i_M)+trans_Comp(71-1,i_S))*cov1[71-1]*ve1[71-1];
				
				n_eff = out_Comp(71,i_V1F);  // number of children effectively protected by MCV1
				n_v1 = out_Comp(71,i_V1S) + out_Comp(71,i_V1I) + out_Comp(71,i_V1R);  // number of children received MCV1
				if (n_v1 > 0.0)
				{
					p_eff = n_eff/n_v1;  // proportion of effective protection among children received MCV1
					if (p_eff > ve2plus) {ve2 = 0.0;}  // vaccine efficacy of 2nd dose conditioned on 1st dose
					else {ve2 = (ve2plus - p_eff)/(1.0 - p_eff);}
				}
				else
				{
					ve2 = ve2plus;
				}
				//Rcout << "\n2nd dose efficacy: " << ve2 << "\n ";
				
				out_Comp(71,i_V2S) = trans_Comp(71,i_V2S)
								- age_w*trans_Comp(71,i_V2S)
								+ age_w*trans_Comp(71-1,i_V2S)*(1.0-cov2)
								+ age_w*trans_Comp(71-1,i_V1S)*cov2*(1.0-ve2);
				
				out_Comp(71,i_V2I) = trans_Comp(71,i_V2I)
								- age_w*trans_Comp(71,i_V2I)
								+ age_w*trans_Comp(71-1,i_V2I)*(1.0-cov2)
								+ age_w*trans_Comp(71-1,i_V1I)*cov2;
				
				out_Comp(71,i_V2R) = trans_Comp(71,i_V2R)
								- age_w*trans_Comp(71,i_V2R)
								+ age_w*trans_Comp(71-1,i_V2R)*(1.0-cov2)
								+ age_w*trans_Comp(71-1,i_V1S)*cov2*ve2
								+ age_w*trans_Comp(71-1,i_V1R)*cov2;
				
				// ve3 = 0, no additional protection for the third dose
				out_Comp(71,i_V3S) = trans_Comp(71,i_V3S)
								- age_w*trans_Comp(71,i_V3S)
								+ age_w*trans_Comp(71-1,i_V3S)
								+ age_w*trans_Comp(71-1,i_V2S)*cov2;
				
				out_Comp(71,i_V3I) = trans_Comp(71,i_V3I)
								- age_w*trans_Comp(71,i_V3I)
								+ age_w*trans_Comp(71-1,i_V3I)
								+ age_w*trans_Comp(71-1,i_V2I)*cov2;
				
				out_Comp(71,i_V3R) = trans_Comp(71,i_V3R)
								- age_w*trans_Comp(71,i_V3R)
								+ age_w*trans_Comp(71-1,i_V3R)
								+ age_w*trans_Comp(71-1,i_V2R)*cov2;
			}
			else 
			{
				out_Comp(a,i_M)   = trans_Comp(a,i_M)
			                      - age_w*trans_Comp(a,i_M)
                			      + age_w*trans_Comp(a-1,i_M)*(1.0-cov1[a-1])
           	    			      - wane*trans_Comp(a,i_M);
           	    
           	    out_Comp(a,i_S)   = trans_Comp(a,i_S)
                                  - age_w*trans_Comp(a,i_S)
                                  + age_w*trans_Comp(a-1,i_S)*(1.0-cov1[a-1])
           	    			      + wane*trans_Comp(a,i_M);
           	    
           	    out_Comp(a,i_I)   = trans_Comp(a,i_I)
                                  - age_w*trans_Comp(a,i_I)
           	    			      + age_w*trans_Comp(a-1,i_I)*(1.0-cov1[a-1]);
           	    
           	    out_Comp(a,i_R)   = trans_Comp(a,i_R)
                                  - age_w*trans_Comp(a,i_R)
                			      + age_w*trans_Comp(a-1,i_R)*(1.0-cov1[a-1]);
           	    
           	    out_Comp(a,i_V1S) = trans_Comp(a,i_V1S)
                                  - age_w*trans_Comp(a,i_V1S)
                	              + age_w*trans_Comp(a-1,i_V1S)
                	              + age_w*(trans_Comp(a-1,i_M)+trans_Comp(a-1,i_S))*cov1[a-1]*(1.0-ve1[a-1]);
           	    
           	    out_Comp(a,i_V1I) = trans_Comp(a,i_V1I)
                                  - age_w*trans_Comp(a,i_V1I)
                	              + age_w*trans_Comp(a-1,i_V1I)
                	              + age_w*trans_Comp(a-1,i_I)*cov1[a-1];
           	    
           	    out_Comp(a,i_V1R) = trans_Comp(a,i_V1R)
                                  - age_w*trans_Comp(a,i_V1R)
                	              + age_w*trans_Comp(a-1,i_V1R)
                	              + age_w*(trans_Comp(a-1,i_M)+trans_Comp(a-1,i_S))*cov1[a-1]*ve1[a-1]
           	    	              + age_w*trans_Comp(a-1,i_R)*cov1[a-1];
           	    
           	    out_Comp(a,i_V1F) = trans_Comp(a,i_V1F)
                                  - age_w*trans_Comp(a,i_V1F)
                                  + age_w*trans_Comp(a-1,i_V1F)
                	              + age_w*(trans_Comp(a-1,i_M)+trans_Comp(a-1,i_S))*cov1[a-1]*ve1[a-1];
           	    
           	    out_Comp(a,i_V2S) = trans_Comp(a,i_V2S) - age_w*trans_Comp(a,i_V2S) + age_w*trans_Comp(a-1,i_V2S);
                out_Comp(a,i_V2I) = trans_Comp(a,i_V2I) - age_w*trans_Comp(a,i_V2I) + age_w*trans_Comp(a-1,i_V2I);
                out_Comp(a,i_V2R) = trans_Comp(a,i_V2R) - age_w*trans_Comp(a,i_V2R) + age_w*trans_Comp(a-1,i_V2R);
                out_Comp(a,i_V3S) = trans_Comp(a,i_V3S) - age_w*trans_Comp(a,i_V3S) + age_w*trans_Comp(a-1,i_V3S);
                out_Comp(a,i_V3I) = trans_Comp(a,i_V3I) - age_w*trans_Comp(a,i_V3I) + age_w*trans_Comp(a-1,i_V3I);
                out_Comp(a,i_V3R) = trans_Comp(a,i_V3R) - age_w*trans_Comp(a,i_V3R) + age_w*trans_Comp(a-1,i_V3R);
            
			//Rcout << a+1 << " ";
			}
		}
		
		// ageing process: 1 week at birth
		pop_fert_SR = 0.0;
		pop_fert_R = 0.0;
		for (int a = 171; a < 189; ++a)  // fertility age: 18-35 years old
		{
			pop_fert_SR += (1.0 - out_Comp(a,i_I) - out_Comp(a,i_V1I) - out_Comp(a,i_V2I) - out_Comp(a,i_V3I))*pop_full[a];
			pop_fert_R += (out_Comp(a,i_R) + out_Comp(a,i_V1R) + out_Comp(a,i_V2R) + out_Comp(a,i_V3R))*pop_full[a];
		}
		prp_R = pop_fert_R/pop_fert_SR;  // Proportion of being born with maternal immunity
		//Rcout << "Proportion of immune = " << prp_R*100 << "%%\n";
	
		out_Comp(0,i_M)   = trans_Comp(0,i_M)   - age_w*trans_Comp(0,i_M)   + age_w*prp_R;
		out_Comp(0,i_S)   = trans_Comp(0,i_S)   - age_w*trans_Comp(0,i_S)   + age_w*(1.0-prp_R);
		out_Comp(0,i_I)   = trans_Comp(0,i_I)   - age_w*trans_Comp(0,i_I)  ;
		out_Comp(0,i_R)   = trans_Comp(0,i_R)   - age_w*trans_Comp(0,i_R)  ;
		out_Comp(0,i_V1S) = trans_Comp(0,i_V1S) - age_w*trans_Comp(0,i_V1S);
		out_Comp(0,i_V1I) = trans_Comp(0,i_V1I) - age_w*trans_Comp(0,i_V1I);
		out_Comp(0,i_V1R) = trans_Comp(0,i_V1R) - age_w*trans_Comp(0,i_V1R);
		out_Comp(0,i_V2S) = trans_Comp(0,i_V2S) - age_w*trans_Comp(0,i_V2S);
		out_Comp(0,i_V2I) = trans_Comp(0,i_V2I) - age_w*trans_Comp(0,i_V2I);
		out_Comp(0,i_V2R) = trans_Comp(0,i_V2R) - age_w*trans_Comp(0,i_V2R);
		out_Comp(0,i_V3S) = trans_Comp(0,i_V3S) - age_w*trans_Comp(0,i_V3S);
		out_Comp(0,i_V3I) = trans_Comp(0,i_V3I) - age_w*trans_Comp(0,i_V3I);
		out_Comp(0,i_V3R) = trans_Comp(0,i_V3R) - age_w*trans_Comp(0,i_V3R);
		out_Comp(0,i_V1F) = trans_Comp(0,i_V1F) - age_w*trans_Comp(0,i_V1F);
		
		
		// =================================================
		// SIA at mid-year
		// =================================================
		// allow up to 5 rounds in a single year
		if ((sia_implement >= 1) && ((t >= (t_start + tstep/2)) && (t < (t_start + tstep/2 + 5)))) 
		{
			IntegerVector alla0 = as<IntegerVector>(siaparm["a0"]);
			IntegerVector alla1 = as<IntegerVector>(siaparm["a1"]);
			NumericVector allsiacov = as<NumericVector>(siaparm["siacov"]);
						
			if ((t - (t_start + tstep/2)) < alla0.size())
			{
				//Rcout << "SIA round: " << ird+1 << "\n";
				int ird = t - (t_start + tstep/2);  // timesteps counting from the beginning of the year
				int a0 = alla0[ird];                // starting target age group for a specific SIA round
				int a1 = alla1[ird];                // ending target age group for a specific SIA round
				double siacov = allsiacov[ird];     // SIA coverage among total population for a specific SIA round
				in_Comp = clone(out_Comp);          // temporary compartments after including transmission and routine vaccination, 254 age groups, 14 states			
					
				double siadose = 0.0;                        // number of total SIA doses
				double pop_0dose = 0.0, pop_non0dose = 0.0;  // number of zero-dose and one-or-more-dose populations 
				
				for (int a = (a0-1); a < a1; ++a) 
				{
					siadose += siacov*pop_full[a];
					pop_0dose += (in_Comp(a,i_M) + in_Comp(a,i_S) + in_Comp(a,i_I) + in_Comp(a,i_R))*pop_full[a];
					pop_non0dose += (1.0 - in_Comp(a,i_M) - in_Comp(a,i_S) - in_Comp(a,i_I) - in_Comp(a,i_R))*pop_full[a];
				}
				//Rcout << "populations: " << pop_0dose << " (zero dose) " << pop_non0dose << " (>= 1 dose)\n";
				if (pop_0dose < 0.5) {pop_0dose = 0.5;}          // avoid zero population in coverage calculation
				if (pop_non0dose < 0.5) {pop_non0dose = 0.5;}

				double siacov1 = 0.0, siacov2 = 0.0; // SIA coverage among zero-dose population and one-or-more-dose population
				
				// SIA implementation based on the weighted logistic function using Portnoy's data
				if (sia_implement == 1)
				{
					double siacov_0dose = exp(-2.621733 + 5.238249*siacov)/(1.0 + exp(-2.621733 + 5.238249*siacov)); // estimated SIA coverage among zero-dose children
					double siadose_0dose = siacov_0dose*pop_0dose;  // number of zero-dose population receiving SIA
					
					if (siadose < siadose_0dose) 
					{
						siacov1 = siadose/pop_0dose;
						siacov2 = 0.0;  
					}
					else 
					{
						if ((siadose - siadose_0dose) >= pop_non0dose)
						{
							siacov1 = (siadose-pop_non0dose)/pop_0dose;
						    siacov2 = 1.0;
						}
						else
						{
							siacov1 = siacov_0dose;
							siacov2 = (siadose - siadose_0dose)/pop_non0dose;
						}
					}
				} 
				else 
				{
				// SIA implementation assuming 7.7% of population are never reached
                    double pr_0dose = pop_0dose/(pop_0dose + pop_non0dose);  // proportion of zero-dose population
					
					if ((pop_non0dose > 0) && ((siacov < 0.923) && (pr_0dose > 0.077)))
					{ // doses given randomly to the population except for the never-reached 
						siacov1 = siadose*((pr_0dose-0.077)/(1-0.077))/pop_0dose;
						siacov2 = siadose*((1-pr_0dose)/(1-0.077))/pop_non0dose;						
					}
					else 
					{ // doses first given to already-vaccinated and then to zero-dose populations				
					    if (siadose > pop_non0dose)
						{
							siacov1 = (siadose-pop_non0dose)/pop_0dose; 
							siacov2 = 1.0;
						} 
						else
						{
							siacov1 = 0.0; 
							siacov2 = siadose/pop_non0dose;
						}
					}
				}
				//Rcout << "SIA coverages: " << siacov1 << " (dose 1) " << siacov2 << " (dose 2)\n";      
				
				if (siacov1 > 1.0) {siacov1 = 1.0;} 
				if (siacov2 > 1.0) {siacov2 = 1.0;} 
				
				
				// campaign vaccination - SIA
				double n_eff = 0.0, n_v1 = 0.0, p_eff = 0.0, ve2 = 0.0;  // initialise for later calculations
				
				for (int a = (a0-1); a < a1; ++a)         // target age groups
				{
					out_Comp(a,i_M)     = in_Comp(a,i_M)
										- in_Comp(a,i_M)*siacov1;     
									
					out_Comp(a,i_S)     = in_Comp(a,i_S)   
										- in_Comp(a,i_S)*siacov1;       
					
					out_Comp(a,i_I)     = in_Comp(a,i_I)   
										- in_Comp(a,i_I)*siacov1; 
									
					out_Comp(a,i_R)     = in_Comp(a,i_R) 
										- in_Comp(a,i_R)*siacov1; 
									
					out_Comp(a,i_V1S)   = in_Comp(a,i_V1S) 
										- in_Comp(a,i_V1S)*siacov2
										+ (in_Comp(a,i_M)+in_Comp(a,i_S))*siacov1*(1.0-ve1[a]);
									
					out_Comp(a,i_V1I)   = in_Comp(a,i_V1I) 
										- in_Comp(a,i_V1I)*siacov2
										+ in_Comp(a,i_I)*siacov1; 
									
					out_Comp(a,i_V1R)   = in_Comp(a,i_V1R) 
										- in_Comp(a,i_V1R)*siacov2
										+ in_Comp(a,i_R)*siacov1
										+ (in_Comp(a,i_M)+in_Comp(a,i_S))*siacov1*ve1[a];
					
					out_Comp(a,i_V1F)   = in_Comp(a,i_V1F) 
										- in_Comp(a,i_V1F)*siacov2
										+ (in_Comp(a,i_M)+in_Comp(a,i_S))*siacov1*ve1[a];	
					
					n_eff = out_Comp(a,i_V1F);  // number of children effectively protected by MCV1
					n_v1 = out_Comp(a,i_V1S) + out_Comp(a,i_V1I) + out_Comp(a,i_V1R);  // number of children received MCV1
					if (n_v1 > 0.0)
					{
						p_eff = n_eff/n_v1;  // proportion of effective protection among children received MCV1
						if (p_eff > ve2plus) {ve2 = 0.0;}  // vaccine efficacy of 2nd dose conditioned on 1st dose
						else {ve2 = (ve2plus - p_eff)/(1.0 - p_eff);}
					}
					else 
					{
						ve2 = ve2plus;
					}
					
					//Rcout << "\n2nd dose efficacy: " << ve2 << "\n ";
					
					out_Comp(a,i_V2S)   = in_Comp(a,i_V2S) 
										- in_Comp(a,i_V2S)*siacov2
										+ in_Comp(a,i_V1S)*siacov2*(1.0-ve2);
								
					out_Comp(a,i_V2I)   = in_Comp(a,i_V2I) 
										- in_Comp(a,i_V2I)*siacov2
										+ in_Comp(a,i_V1I)*siacov2;
								
					out_Comp(a,i_V2R)   = in_Comp(a,i_V2R) 
										- in_Comp(a,i_V2R)*siacov2
										+ in_Comp(a,i_V1S)*siacov2*ve2
										+ in_Comp(a,i_V1R)*siacov2;							
					
					// ve3 = 0, no additional protection for the third dose		
					out_Comp(a,i_V3S)   = in_Comp(a,i_V3S) 
										+ in_Comp(a,i_V2S)*siacov2;	
								
					out_Comp(a,i_V3I)   = in_Comp(a,i_V3I) 
										+ in_Comp(a,i_V2I)*siacov2;
								
					out_Comp(a,i_V3R)   = in_Comp(a,i_V3R) 
										+ in_Comp(a,i_V2R)*siacov2;
					
					//Rcout << a+1 << " ";
				}
			}
		}
					
	    in_Comp = clone(out_Comp);
    }

	outp["cases"] = newinfect;
    outp["out_Comp"] = in_Comp;

  return outp;
}