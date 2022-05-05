#include "fitter.h"
#include <sstream>
#include <iostream>
#include <cmath>




using namespace std;

vector<vector<double>>  Read::prin_corr(string file_location ){
        	
  vector<vector<double>> prin;
  vector<double> prin_tmp;
  ifstream inFile; 
  inFile.open(file_location);
  int Ncfgs=0; int Nt = 0;
  if(inFile.is_open()){
	
	  
	  string line;
	  getline(inFile,line);
	  istringstream iss(line);
	  int cfg_tmp; int nt_tmp;
	  iss >> cfg_tmp >> nt_tmp;
	  Ncfgs+=cfg_tmp;Nt+=nt_tmp;

	
	  while(getline(inFile,line)){
		
	    istringstream iss(line);
            int t;double val;
	    iss >> t >> val;
			
	    prin_tmp.push_back( val);
						
			} // while
	inFile.close();} // file open

	else{cerr<<"can't find input file: "<< file_location<<endl;

            }
       cout<<"Reshaping correlator "<<prin_tmp.size()<<endl;
       
       
       for( int t = 0; t<Nt; t++){
       vector<double> vec;
         for(int cfg=0; cfg<Ncfgs;cfg++){

	   vec.push_back(prin_tmp[Nt*cfg+t]);
	                       } // cfg
	prin.push_back(vec);
                                          } // t
        
       
       return prin;

};

itpp::mat pseudoinvertSVD(const itpp::mat& A, double tol)
{
  int M = A.rows();
  int N = A.cols();

  itpp::mat U(M,M);
  itpp::mat V(N,N);
  itpp::vec s(M);

  itpp::svd(A, U, s, V);

  //  cout << "done the svd" << endl;

  double smax = itpp::max(s);
  int count = 0;

  itpp::mat Sinv(N,M);
  Sinv.zeros();
  for(int i = 0; i < N; i++)
    { //test for small singular values - kills those smaller than tol of maximum
      if(s(i) > tol * smax){ Sinv(i,i) = 1./s(i);}
      else{Sinv(i,i) = 0; count++;}
    }

  //nReset = count;
   cout << "reset " << count << " sing values" << endl;
  itpp::mat inverse = V*Sinv*transpose(U);
  return inverse;
};


vector<double> Ensemble::average(vector<vector<double>> prin_corr){

	vector<double> prin_corr_avg;
        int Nt = prin_corr.size(); int Ncfgs = prin_corr[0].size();
	for(int t =0; t < Nt; t++){
		double avg = 0;
		for( int cfg = 0; cfg<Ncfgs;cfg++){
			avg+=prin_corr[t][cfg];

		} // cfg
		
		prin_corr_avg.push_back(avg/Ncfgs);

	} //t 
	return prin_corr_avg;

};
vector<double> Ensemble::peek(vector<vector<double>> prin_corr, int cfg){
	vector<double> prin_corr_cfg;
	int Nt = prin_corr.size();
	for (int t=0; t<Nt;t++){prin_corr_cfg.push_back(prin_corr[t][cfg]);}	
	return prin_corr_cfg;

}

void Ensemble::scale_up(vector<vector<double>> prin_corr){

        
	int Ncfgs = prin_corr[0].size();
	
	if (Ncfgs < 2){cerr<<"Need more than one configuration"<<endl;}
	
	int Nt = prin_corr.size();

	vector<double>	prin_means =  average (prin_corr) ;

	for( int t = 0; t<Nt; t++){
		for( int cfg = 0; cfg<Ncfgs; cfg++){
			prin_corr[t][cfg] = (Ncfgs /( Ncfgs-1.0) )*prin_means[t] - (1.0/(Ncfgs-1.0))*prin_corr[t][cfg];

		}//cfg

	}//t
	
};

itpp::mat Ensemble::invcovariance(int tmin, int tmax, double tol, vector<vector<double>> prin_corr){

	int Nt = prin_corr.size();
	int Ncfgs = prin_corr[0].size();
        cout<<"Number of configs "<<Ncfgs<<endl;
        itpp::mat cova;

	cova.set_size(tmax-tmin, tmax-tmin);
        
	for( int t1=0; t1<tmax-tmin;t1++){
		vector<double> tmp;
		for( int t2 = t1; t2< tmax-tmin; t2++){
		  
		  double cov = 0;
		  for(int cfg =0;cfg<Ncfgs;cfg++){


		  	  cov+=(prin_corr[t1+tmin][cfg]-average(prin_corr)[t1+tmin])*(prin_corr[t2+tmin][cfg]-average(prin_corr)[t2+tmin]);
		  } // cfg
		cova(t1,t2) = cov/Ncfgs;
		cova(t2,t1) = cov/Ncfgs;
		} //t2
	  
	}//t1
        itpp::mat cinv = pseudoinvertSVD(cova, tol);
	return cinv;
};


vector<double> Fit_stuff::double_exp(int t0, int tmin, int tmax, double pars[3]){

		double A = pars[0]; double m0= pars[1]; double m1= pars[2];

		vector<double> lambda;

		for( int t=tmin; t<tmax;t++){

			lambda.push_back((1.0-A)*exp(-m0*(t-t0))+A*exp(-m1*(t-t0)));
		}
		return lambda;

};
vector<double> Fit_stuff::double_exp_jac(int t0, int t, double pars[3]){


		double A = pars[0]; double m0 = pars[1]; double m1 = pars[2];
	
		vector<double> jac;
		double jac_A = - exp(-m0 * (t-t0) ) + exp( -m1 * (t-t0) );
		double jac_m0 = (1.0-A) * (-t + t0)* exp(-m0*(t-t0));
		double jac_m1 = -A * (t - t0)* exp(-m1*(t-t0));
		jac.push_back(jac_A);jac.push_back(jac_m0);jac.push_back(jac_m1);
		return jac;

};

double chisquare(const gsl_vector *v, void * params){
        double A, m0, m1;
	A = gsl_vector_get(v,0); m0 = gsl_vector_get(v,1); m1 = gsl_vector_get(v,2);
	corr_pack *corr = (corr_pack *) params;
	itpp:: mat invcov = corr->invcov;
        vector<double> prin_corr = corr->prin_corr;
	double total = 0;
	int t0 = corr->t0;
	int tmin = corr->tmin;
	int tmax = corr->tmax;
        
	for (int t1=0;t1<tmax-tmin;t1++){
		for ( int t2=0; t2<tmax-tmin;t2++){
		total+=(( (1.0-A)*exp(-m0*(t1+tmin-t0))+A*exp(-m1*(t1+tmin-t0))) -prin_corr[t1+tmin] ) * invcov(t1,t2)* (( (1.0-A)*exp(-m0*(t2+tmin-t0))+A*exp(-m1*(t2+tmin-t0))) -prin_corr[t2+tmin] );
	}}
	cout<<total<<endl;
	return total/(tmax-tmin-3.0);
};

void chisquare_grad(const gsl_vector *v, void * params, gsl_vector *df){
        double A, m0, m1;
	Fit_stuff fit;
        A = gsl_vector_get(v,0); m0 = gsl_vector_get(v,1); m1 = gsl_vector_get(v,2);
        corr_pack *corr =(corr_pack *) params;
        itpp:: mat invcov = corr->invcov;
        vector<double> prin_corr = corr->prin_corr;
        double total = 0;
        int t0 = corr->t0;
        int tmin = corr->tmin;
        int tmax = corr->tmax;
        double pars_new[3];
	pars_new[0] = A; pars_new[1] = m0; pars_new[2] = m1;

	for(int i = 0; i < 3;i++){
         
          double total = 0;
          vector<double> lambda = fit.double_exp(t0, tmin, tmax, pars_new);
          for (int t1=0; t1<tmax-tmin;t1++){vector<double> jac1 = fit.double_exp_jac(t0, t1+tmin, pars_new);
          for(int t2 =0; t2<tmax-tmin;t2++){
                vector<double> jac2 = fit.double_exp_jac(t0, t2+tmin, pars_new);


                total+=( (jac1[i])*invcov(t1,t2)*(lambda[t2] - prin_corr[tmin+t2]) + (lambda[t1] - prin_corr[tmin+t1])*invcov(t1,t2)*(jac2[i]))/(tmax-tmin-3.0) ;
                                }//t2

        }//t1
          
	gsl_vector_set(df,i,total);
        }

};

void chisquare_fdf(const gsl_vector *x, void * params, double *f, gsl_vector *df){

*f = chisquare(x,params);
chisquare_grad(x,params,df);
}
/*double Fit_stuff::chisq(vector<vector<double>> prin_corr, itpp::mat invcov, int cfg, const params& pars, const fit_controls& controls ){

	int Nt = prin_corr.size();
	int tmin = controls.tmin;
	int tmax = controls.tmax;
	int t0 = controls.t0;
	double tol = controls.svd_cutoff;
	double chisq = 0;
	vector<double> lambda = double_exp(t0, tmin, tmax, pars);
	Ensemble ensem;
        for (int t1=0; t1<tmax-tmin;t1++){
	  for(int t2 =0; t2<tmax-tmin;t2++){
	       chisq+=	(lambda[t1] - prin_corr[tmin+t1][cfg])*invcov(t1,t2)*(lambda[t2] - prin_corr[tmin+t2][cfg]);
	}// t2
	}//t1
	cout << chisq << endl;
	return chisq/(tmax-tmin-3.0) ;
};*/

params Fit_stuff::cg_minizer(vector<vector<double>> prin_corr, int cfg, params& pars, const fit_controls& controls ){

        int Nt = prin_corr.size();
        int tmin = controls.tmin;
        int tmax = controls.tmax;
        int t0 = controls.t0;
        double tol = controls.svd_cutoff;
//	int iter = controls.iterations;
        double eta = controls.eta;
//	params pars_new = pars;

        
	Ensemble ensem;
        cout<<"Making covariance"<<endl;
	itpp::mat invcov = ensem.invcovariance(tmin, tmax, tol, prin_corr);
	// get derivative of lambdas
	cout<<"Beginning fitting "<<endl;
	
	
	         

         size_t iter = 0;
	 int status;
 
	 const gsl_multimin_fdfminimizer_type *T;
	 gsl_multimin_fdfminimizer *s;
         //void * CORR;
	 corr_pack  CORR;

         
	 CORR.invcov = invcov;
	
	 CORR.t0 = t0; CORR.tmin = tmin; CORR.tmax = tmax;
	 
         CORR.prin_corr = ensem.peek(prin_corr,cfg);
	 corr_pack * corr_ = &CORR;
	 gsl_vector * x;
	 gsl_multimin_function_fdf minimize;
	 minimize.n = 3;
	 minimize.f = &chisquare;
	 minimize.df = &chisquare_grad;
	 minimize.fdf = &chisquare_fdf;
	 minimize.params = corr_;
	// minimize.params
	 //corr_pack *corr = (corr_pack *)params;
	
	 x = gsl_vector_alloc(3);
	 gsl_vector_set(x,0,0.2); gsl_vector_set(x,1,0.77); gsl_vector_set(x,2,1.0);

	 T = gsl_multimin_fdfminimizer_conjugate_fr;
	 s = gsl_multimin_fdfminimizer_alloc(T,3);
	 gsl_multimin_fdfminimizer_set(s, &minimize, x, 0.01,1.0e-4);
	 do
	 {
		 iter++;
		 status = gsl_multimin_fdfminimizer_iterate(s);
		 if(status)
		    break;
		 
			status = gsl_multimin_test_gradient(s->gradient,1.e-3);

		if(status == GSL_SUCCESS){printf("Minimum found at:\n");}
                  
              cout<< iter << " "<< gsl_vector_get (s->x, 0) <<" " << gsl_vector_get (s->x, 1)<<" "<<gsl_vector_get(s->x, 2) << endl;
	 }
	         while (status == GSL_CONTINUE && iter < 1.e4);

		 gsl_multimin_fdfminimizer_free(s);
		 gsl_vector_free(x);
	 
	 /*	vector<double> grad;
	for(int i = 0; i < 3;i++){	
	  
	  double total = 0;
	  vector<double> lambda = double_exp(t0, tmin, tmax, pars_new);
	  for (int t1=0; t1<tmax-tmin;t1++){vector<double> jac1 = double_exp_jac(t0, t1+tmin, pars_new);
          for(int t2 =0; t2<tmax-tmin;t2++){
		vector<double> jac2 = double_exp_jac(t0, t2+tmin, pars_new); 
                
            
		total+=( (jac1[i])*invcov(t1,t2)*(lambda[t2] - prin_corr[tmin+t2][cfg]) + (lambda[t1] - prin_corr[tmin+t1][cfg])*invcov(t1,t2)*(jac2[i]))/(tmax-tmin-3.0) ;
			       	}//t2

	}//t1 
	  grad.push_back(total);

	} // grad
        double chi  = chisq(prin_corr, invcov, cfg, pars_new, controls);
	//pars_new.A = pars_new.A - eta* chi / grad[0]; 
	pars_new.m0 = pars_new.m0 - eta* chi / grad[1];
       //	pars_new.m1 = pars_new.m1 - eta*chi/grad[2];

	cout<<"A= "<<pars_new.A<<" m0= "<<pars_new.m0<<" m1= "<<pars_new.m1<<endl;
	*/


//	return pars_new;
	};



int main(){
        
	
	Read rd;
        vector<vector<double>> corr = rd.prin_corr("prin_correlator_1.jack");
        Ensemble ensem;

	ensem.scale_up(corr);

        Fit_stuff fit;
	fit_controls controls;
        controls.t0 = 3;
        controls.tmin =5 ;
        controls.tmax = 15;
        controls.svd_cutoff = 1.e-6;
        controls.iterations = 10000;
        controls.eta = 1.e-1;

        int cfg = 1;
	params pars;
	pars.A = 0.2; pars.m0 = 0.77; pars.m1 = 1.0;
        fit.cg_minizer(corr, cfg, pars, controls);
return 0;
}
