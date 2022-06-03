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

itpp::mat pseudoinvertSVD(const itpp::mat& A, double tol, int& nReset)
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

  nReset = count;
  // cout << "reset " << count << " sing values" << endl;
  itpp::mat inverse = V*Sinv*transpose(U);
  return inverse;
};


itpp::mat invertSVDNorm(const itpp::mat& A, double tol, int& nReset)
{
  //normalise the matrix to O(1) using the sqrt of the diag
  //only for square matrices now
  if(A.rows() != A.cols()){cerr << "square matrices only in invertSVDNorm" << endl; exit(1);}

  itpp::mat D(A.rows(), A.cols());
  D.zeros();
  for(int i=0; i < A.rows(); i++){D(i,i) = 1.0 / sqrt(fabs(A(i,i)));}

  itpp::mat DAD = D*A*D;

  // cout << "normalised the covariance" << endl;

  itpp::mat AA = pseudoinvertSVD(DAD, tol, nReset);

  return D*AA*D;
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
double Ensemble::average(vector<double> mass){


	double m_jack = 0;
	for(int cfg =0;cfg< mass.size();cfg++){m_jack+=mass[cfg];}
	return m_jack/mass.size();


};

double Ensemble::variance(vector<double> mass){

        int Ncfgs = mass.size();
        double avg = average(mass);
        double var =0;
        for(int cfg =0;cfg<Ncfgs;cfg++){

                var+=pow(mass[cfg]-avg,2);

        }
        return var/Ncfgs;

};

vector<double> Ensemble::scale_down(vector<double> mass){

	vector<double> m_jack;
	int Ncfgs = mass.size();

        if (Ncfgs < 2){cerr<<"Need more than one configuration"<<endl;}

        

        double av = average(mass);

      
                for( int cfg = 0; cfg<Ncfgs; cfg++){
                        double jack = Ncfgs*av - (Ncfgs-1.0)*mass[cfg];
                        m_jack.push_back(jack);
                }//cfg
                

     
        return m_jack;

};

vector<double> Ensemble::peek(vector<vector<double>> prin_corr, int cfg){
	vector<double> prin_corr_cfg;
	int Nt = prin_corr.size();
	for (int t=0; t<Nt;t++){prin_corr_cfg.push_back(prin_corr[t][cfg]);}	
	return prin_corr_cfg;

}

vector<vector<double>> Ensemble::scale_up(vector<vector<double>> prin_corr){

        vector<vector<double>> p_corr; 
	int Ncfgs = prin_corr[0].size();
	
	if (Ncfgs < 2){cerr<<"Need more than one configuration"<<endl;}
	
	int Nt = prin_corr.size();

	vector<double>	prin_means =  average (prin_corr) ;

	for( int t = 0; t<Nt; t++){vector<double> tmp;
		for( int cfg = 0; cfg<Ncfgs; cfg++){
			double jack = (Ncfgs /( Ncfgs-1.0) )*prin_means[t] - (1.0/(Ncfgs-1.0))*prin_corr[t][cfg];
			tmp.push_back(jack);
		}//cfg 
		p_corr.push_back(tmp);

	}//t
	return p_corr;
	
};

vector<vector<double>> Ensemble::scale_down(vector<vector<double>> prin_corr){

        vector<vector<double>> p_corr;
        int Ncfgs = prin_corr[0].size();

        if (Ncfgs < 2){cerr<<"Need more than one configuration"<<endl;}

        int Nt = prin_corr.size();

        vector<double>  prin_means =  average (prin_corr) ;

        for( int t = 0; t<Nt; t++){vector<double> tmp;
                for( int cfg = 0; cfg<Ncfgs; cfg++){
                        double jack = Ncfgs*prin_means[t] - (Ncfgs-1.0)*prin_corr[t][cfg];
                        tmp.push_back(jack);
                }//cfg 
                p_corr.push_back(tmp);

        }//t
        return p_corr;

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
		  cov/=(Ncfgs);//*(Ncfgs-1.0));
		  cov/=(Ncfgs-1.0);
  		  cova(t1,t2) = cov;

		cova(t2,t1) = cov;
		} //t2
	  
	}//t1
	int nResetCovSingVals;
        itpp::mat cinv = invertSVDNorm(cova, tol, nResetCovSingVals); 
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
Fit_stuff fit;
Ensemble ensem;
double Fit_stuff::get_chisq(vector<vector<double>> prin_corr, vector<double> lambda, itpp::mat invcov, int tmin, int tmin_fit, int tmax){
	
	int Ncfgs = prin_corr[0].size();
	double total = 0.;
	vector<double> p_mean = ensem.average(prin_corr);	
	for (int t1=0;t1<tmax-tmin_fit;t1++){
                for ( int t2=0; t2<tmax-tmin_fit;t2++){

                total+=(lambda[t1] -p_mean[t1+tmin_fit] ) *invcov(t1+tmin_fit-tmin,t2+tmin_fit-tmin)* (lambda[t2]-p_mean[t2+tmin_fit] );

		//cout << t1<<" "<<t2<<" :"<<(lambda[t1] -p_mean[t1+tmin_fit] ) <<" "<< invcov(t1+tmin_fit-tmin,t2+tmin_fit-tmin) <<" "<< (lambda[t2]-p_mean[t2+tmin_fit] ) << endl;
		}}
	



        return total/(tmax-tmin_fit-3);

	
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
	int tmin_fit = corr->tmin_fit;
	

	double pars_new[3];
        pars_new[0] = A; pars_new[1] = m0; pars_new[2] = m1;


        vector<double> lambda = fit.double_exp(t0, tmin_fit, tmax, pars_new);
	for (int t1=0;t1<tmax-tmin_fit;t1++){
		for ( int t2=0; t2<tmax-tmin_fit;t2++){
         
		total+=(lambda[t1] -prin_corr[t1+tmin_fit] ) *invcov(t1+tmin_fit-tmin,t2+tmin_fit-tmin)* (lambda[t2]-prin_corr[t2+tmin_fit] );
	}}
	
	return total/(tmax-tmin_fit-3.0);
};

void chisquare_grad(const gsl_vector *v, void * params, gsl_vector *df){
        double A, m0, m1;

        A = gsl_vector_get(v,0); m0 = gsl_vector_get(v,1); m1 = gsl_vector_get(v,2);
        corr_pack *corr =(corr_pack *) params;
        itpp:: mat invcov = corr->invcov;
        vector<double> prin_corr = corr->prin_corr;
        double total = 0;
        int t0 = corr->t0;
        int tmin = corr->tmin;
        int tmax = corr->tmax;
	int tmin_fit = corr-> tmin_fit;
        double pars_new[3];
	pars_new[0] = A; pars_new[1] = m0; pars_new[2] = m1;

	for(int i = 0; i < 3;i++){
         
          double total = 0;
          vector<double> lambda = fit.double_exp(t0, tmin_fit, tmax, pars_new);
          for (int t1=0; t1<tmax-tmin_fit;t1++){vector<double> jac1 = fit.double_exp_jac(t0, t1+tmin_fit, pars_new);
		  
          for(int t2 =0; t2<tmax-tmin_fit;t2++){
                vector<double> jac2 = fit.double_exp_jac(t0, t2+tmin_fit, pars_new);
	//	for(int cfg=0;cfg<prin_corr[0].size();cfg++){

                total+=( (jac1[i])*invcov(t1+tmin_fit-tmin,t2+tmin_fit-tmin)*(lambda[t2] - prin_corr[tmin_fit+t2]) + (lambda[t1] - prin_corr[tmin_fit+t1])*invcov(t1+tmin_fit-tmin,t2+tmin_fit-tmin)*(jac2[i]))/(tmax-tmin_fit-3.0) ;
	//	} // cfg                
		}//t2

        }//t1
          
	gsl_vector_set(df,i,total);
	
        }

};

void chisquare_fdf(const gsl_vector *x, void * params, double *f, gsl_vector *df){

*f = chisquare(x,params);
chisquare_grad(x,params,df);
}




vector<double> Fit_stuff::cg_minizer(vector<vector<double>> prin_corr,itpp::mat invcov,int cfg,  int t0, int tmin, int tmax, int tmin_fit ){
	 
         int Nt = prin_corr.size();
     
         vector<double> output;
        	         
         size_t iter = 0;
	 int status;
 
	 const gsl_multimin_fdfminimizer_type *T;
	 gsl_multimin_fdfminimizer *s;
         //void * CORR;
	 corr_pack  CORR;
         
	 CORR.invcov = invcov;
	
	 CORR.t0 = t0; CORR.tmin = tmin; CORR.tmax = tmax; CORR.tmin_fit = tmin_fit;
	 
         CORR.prin_corr = ensem.peek(prin_corr,cfg);
	// CORR.prin_corr = prin_corr;
	 corr_pack * corr_ = &CORR;
	 gsl_vector * x;
	 gsl_multimin_function_fdf minimize;
	 minimize.n = 3;
	 minimize.f = &chisquare;
	 minimize.df = &chisquare_grad;
	 minimize.fdf = &chisquare_fdf;
	 minimize.params = corr_;


	 x = gsl_vector_alloc(3);
	 gsl_vector_set(x,0,0.506); gsl_vector_set(x,1,0.5009); gsl_vector_set(x,2,1.7547);

	 T = gsl_multimin_fdfminimizer_conjugate_fr;
	 s = gsl_multimin_fdfminimizer_alloc(T,3);
	 gsl_multimin_fdfminimizer_set(s, &minimize, x, 1.e-1,1.0e-3);
	 do
	 {
		 iter++;
		 status = gsl_multimin_fdfminimizer_iterate(s);
		 if(status)
		    break;
		 
			status = gsl_multimin_test_gradient(s->gradient,1.e-1);

		if(status == GSL_SUCCESS){
			//cout<<"Minimum found "<<iter<<endl;
			output.push_back(gsl_vector_get (s->x, 0)); output.push_back(gsl_vector_get (s->x, 1));output.push_back(gsl_vector_get (s->x, 2));}
                  
             // cout<< iter << " "<< gsl_vector_get (s->x, 0) <<" " << gsl_vector_get (s->x, 1)<<" "<<gsl_vector_get(s->x, 2) << endl;
	 }
	         while (status == GSL_CONTINUE && iter < 1.e4);
	         output.push_back(gsl_vector_get (s->x, 0)); output.push_back(gsl_vector_get (s->x, 1));output.push_back(gsl_vector_get (s->x, 2));
		 gsl_multimin_fdfminimizer_free(s);
		 gsl_vector_free(x);
	         return output;


	};

fit_params Fit_stuff::fit_prin_corr(vector<vector<double>> p_corr_jack, itpp::mat invcov, int t0, int tmin, int tmax, int tslice_num){
	
        //vector<vector<double>> p_corr_jack = ensem.scale_up(corr);
	int Ncfgs = p_corr_jack[0].size();
        vector<double> A;vector<double> m0;vector<double>m1;


        for(int cfg =0;cfg<Ncfgs;cfg++){
        	
		vector<double> answer = cg_minizer(p_corr_jack, invcov, cfg,  t0, tmin, tmax, tslice_num);
                A.push_back(answer[0]);m0.push_back(answer[1]);m1.push_back(answer[2]);
	}
	 
	vector<double> A_jack = ensem.scale_down(A); vector<double> m0_jack = ensem.scale_down(m0); vector<double> m1_jack = ensem.scale_down(m1);

	fit_params result;
	result.A = A_jack; result.m0 = m0_jack; result.m1 = m1_jack;
	// get chisq
	
	double pars_new[3];
        pars_new[0] = ensem.average(A_jack); pars_new[1] = ensem.average(m0_jack); pars_new[2] = ensem.average(m1_jack);


        vector<double> lambda = fit.double_exp(t0, tslice_num, tmax, pars_new);
	vector<vector<double>> p_corr = ensem.scale_down(p_corr_jack);


	result.chisq = get_chisq(p_corr_jack, lambda, invcov, tmin, tslice_num, tmax);
	return result;
	
};




 int main(int argc, char *argv[]){

  if(argc != 7){ cerr << "usage: <filename> <t0> <tmin> <tmax> <min_tslices> <svdcut>" << endl; exit(1); }
  string filen; {istringstream val(argv[1]); val >> filen;}
  int t0; {istringstream val(argv[2]); val >> t0;}
  int tmin; {istringstream val(argv[3]); val >> tmin;}
  int tmax; {istringstream val(argv[4]); val >> tmax;}
  int minTSlices; {istringstream val(argv[5]); val >> minTSlices;}
  double svd_cutoff; {istringstream val(argv[6]); val >> svd_cutoff;}       
	
	Read rd;
        vector<vector<double>> corr = rd.prin_corr(filen);
        

//	ensem.scale_up(corr);

        Fit_stuff fit;


        int Ncfgs = corr[0].size();
        
	
	itpp::mat invcov = ensem.invcovariance(tmin, tmax, svd_cutoff, corr);
 
	vector<vector<double>> p_corr_jack = ensem.scale_up(corr);

#ifdef _OPENMP
  int nthr = omp_get_max_threads();
  
  
#pragma omp parallel for num_threads(nthr)  default(shared)

	for(int t=minTSlices; t<tmax-tmin;t++){
	fit_params result = fit.fit_prin_corr(p_corr_jack, invcov, t0, tmin, tmax, t);
	//cout<<"t: "<<t<<"    "<<omp_get_thread_num()  <<endl;
	cout<<ensem.average(result.m0)<<" "<<sqrt(ensem.variance(result.m0)/(Ncfgs-1))<<"    CHISQ= "<<result.chisq<<endl;
	}
	
	
#endif
	return 0;
}
