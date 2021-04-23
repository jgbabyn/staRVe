// Headers are included in staRVe.cpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

//Little helper for time 
template<class Type>
Type rhoTrans(Type x){
  return Type(2)/(Type(1) + exp(-Type(2)*x))-Type(1);
}

//For barrier model SPDE
template<class Type>
Type mayo_on_sundae(objective_function<Type>* obj) {

  //Call in my friends to make this work
  using namespace R_inla_barrier;
  using namespace density;
  
  // Read in data / parameters / random effectos from R
  DATA_INTEGER(n_time); //not really needed?
  DATA_INTEGER(distribution_code);
  DATA_INTEGER(link_code);

  DATA_VECTOR(obs_y);
  DATA_VECTOR_INDICATOR(keep,obs_y);
  DATA_MATRIX(mean_design);
  DATA_IVECTOR(sample_size);
  
  //Get da FEM matrices
  DATA_STRUCT(fem,fem_barrier_t);
  //The Projector Matrix A for points not on nodes, if any
  DATA_SPARSE_MATRIX(A);

  //Parameters
  PARAMETER_VECTOR(working_response_pars); //Response dist parameters except for mean
  PARAMETER_VECTOR(mean_pars); // Fixed effects, just Beta...
  PARAMETER(log_range); //range for Matern of SPDE, barrier one is just fraction of this 
  PARAMETER(log_sigma_u); //marginal variance of Matern
  PARAMETER_ARRAY(Z); //holds random effects from SPDE
  PARAMETER_VECTOR(logit_time_ar1); // AR1 time rho thing, for now unused
  PARAMETER_VECTOR(log_time_sd); //AR1 scale factor?
  

  //Some transforms
  vector<Type> ranges(2);
  ranges(0) = exp(log_range);
  ranges(1) = ranges(0)*0.1;
  Type sigma_u = exp(log_sigma_u);
  Type time_sd;
  Type rho;

  //If there is time stuff
  if(log_time_sd.size() != 0){
    time_sd = exp(log_time_sd(0));
  }
  if(logit_time_ar1.size() != 0){
    rho = rhoTrans(logit_time_ar1(0));
  }
    // Convert parameters from working scale to natural scale

  // Cross-check distribution_code order with family.hpp
  //You should probably shove this in a function now that I use it again
  vector<Type> response_pars = working_response_pars;
  switch(distribution_code) {
    case 0 : response_pars(0) = exp(working_response_pars(0)); break; // Normal, sd>0
    case 1 : break; // Poisson, NA
    case 2 : response_pars(0) = exp(working_response_pars(0))+1; break; // Neg. Binom., overdispersion > 1
    case 3 : break; // Bernoulli, NA
    case 4 : response_pars(0) = exp(working_response_pars(0)); break; // Gamma, sd>0
    case 5 : response_pars(0) = exp(working_response_pars(0)); break; // Log-Normal, sd>0
    case 6 : break; // Binomial, NA
    case 7 : break; // AtLeastOneBinomai, NA
    case 8 : response_pars(0) = exp(working_response_pars(0)); break; // Conway-Maxwell-Poisson, dispersion > 0
    default : response_pars(0) = exp(-1*working_response_pars(0)); break; // Normal, sd>0
  }

  // Set up response dist and link function
  inv_link_function inv_link = {link_code}; //Why brackets?
  response_density y_density = {distribution_code};
  glm<Type> family(inv_link,
		   y_density,
		   mean_pars,
		   response_pars);
  
  Type nll = 0.0;

  //Make the precision matrix Q
  SparseMatrix<Type> Q = Q_barrier(fem,ranges,sigma_u);

  //Timey wimey bits
  GMRF_t<Type> spat(Q);

  if(log_time_sd.size() != 0){
    N01<Type> nllN01;
    AR1_t<N01<Type> > time(rho,nllN01);
    SCALE_t<AR1_t<N01<Type> > > tscaled(time,time_sd);
    SEPARABLE_t<SCALE_t<AR1_t<N01<Type> > >, GMRF_t<Type> > spattime(tscaled,spat);
    nll += spattime(Z);
  }else{ //Just space
    nll += spat(Z);
  }
  
  
  //add it's contribution to the nll
  

  matrix<Type> AZ(A.rows(),Z.cols());
  for(int i = 0; i < Z.cols();i++){
    AZ.col(i) = A*Z.matrix().col(i);
  }

  vector<Type> eta(obs_y.size());
  for(int i = 0; i < mean_design.rows();i++){
    eta(i) = family.inv_link(mean_design.row(i),AZ(i)); //has to be column order I guess
   }

  for(int i = 0; i < obs_y.size(); i++){
     nll -= keep(i)*family.log_density(obs_y(i),eta(i), sample_size(i));
   }
  
  
  
  return nll;
}

//For regular SPDE?
template<class Type>
Type spde_classic_model(objective_function<Type>* obj) {

  Type nll = 0.0;

  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
