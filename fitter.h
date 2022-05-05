#ifndef __FITTER_H__
#define __FITTER_H__

#include <vector>
#include <fstream>
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

struct corr_pack{

	vector<double> prin_corr; 
	itpp::mat invcov; 
	int t0;
        int tmin;
        int tmax;
};

class Ensemble{

    public: 
            vector<double> average( vector<vector<double>> prin_corr);
	    
	    void scale_up ( vector<vector<double>> prin_corr );
            vector<double> peek(vector<vector<double>> prin_corr, int cfg);  
	    itpp::mat invcovariance(int tmin, int tmax, double tol, vector<vector<double>> prin_corr); 
};

class Read{

    public:
	    vector<vector<double>> prin_corr( string file_location ); // read and return principle correlation function


};


class Fit_stuff{

	
	public:
		vector<double> double_exp(int t0, int tmin, int tmax, double pars[3]);
		vector<double>	double_exp_jac(int t0, int t, double pars[3]);
                double chisq(vector<vector<double>> prin_corr, itpp::mat invcov, int cfg, const params& pars, const fit_controls& controls );
		params cg_minizer(vector<vector<double>> prin_corr, int cfg, params& pars, const fit_controls& controls );

};




#endif
