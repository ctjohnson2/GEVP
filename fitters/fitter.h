#ifndef __FITTER_H__
#define __FITTER_H__

#include <vector>
#include <fstream>
#include <omp.h>
#include <iostream>
#include <itpp/stat/misc_stat.h>
#include <itpp/base/timing.h>
#include <itpp/base/mat.h>
#include <itpp/itbase.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_vector.h>

using namespace std;

struct params{

                double A; double m0; double m1;
        };

        struct fit_controls{
        int t0;
        int tmin;
        int tmax;
        double svd_cutoff;
        int iterations;
        double eta;
        };

struct fit_params{

	vector<double> A; vector<double> m0; vector<double> m1; double chisq;
};
struct corr_pack{

	vector<double> prin_corr; 
	itpp::mat invcov; 
	int t0;
        int tmin;
        int tmax;
	int tmin_fit;
};

class Ensemble{

    public: 
            vector<double> average( vector<vector<double>> prin_corr);
	    double average(vector<double> mass);
	    double variance(vector<double> mass);
	    vector<vector<double>> scale_up ( vector<vector<double>> prin_corr );
	    vector<double> scale_down(vector<double> mass);
            vector<vector<double>> scale_down ( vector<vector<double>> prin_corr );
	    vector<double> peek(vector<vector<double>> prin_corr, int cfg);  
	    itpp::mat invcovariance(int tmin, int tmax, double tol, vector<vector<double>> prin_corr); 
};

class Read{

    public:
	    vector<vector<double>> prin_corr( string file_location ); // read and return principle correlation function
            void write_mass_jack(vector<double> mass_jack, string name, int Ncfgs);

};


class Fit_stuff{

	
	public:
		vector<double> double_exp(int t0, int tmin, int tmax, double pars[3]);
		vector<double>	double_exp_jac(int t0, int t, double pars[3]);
                double chisq(vector<vector<double>> prin_corr, itpp::mat invcov, int cfg, const params& pars, const fit_controls& controls );
		vector<double> cg_minizer(vector<vector<double>> prin_corr,itpp::mat invcov,int cfg,  int t0, int tmin, int tmax, int tslice_num,double init[3] );
		fit_params fit_prin_corr(vector<vector<double>> p_corr_jack, itpp::mat invcov, int t0, int tmin, int tmax, int tslice_num);
		double get_chisq(vector<vector<double>> prin_corr, vector<double> lambda, itpp::mat invcov, int tmin, int tmin_fit, int tmax);
};




#endif
