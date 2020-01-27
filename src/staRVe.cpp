#include <TMB.hpp>
#include <iostream>
using namespace density;

#include "include/data_in.hpp"
#include "include/time_segment.hpp"
#include "include/family.hpp"
#include "include/glm.hpp"

#include "include/covariance.hpp"
#include "include/kriging.hpp"
#include "include/nngp.hpp"
#include "include/observations.hpp"




template<class Type>
Type objective_function<Type>::operator() () {
  DATA_INTEGER(n_time);
  DATA_INTEGER(distribution_code);
  DATA_INTEGER(link_code);

  DATA_IVECTOR(y_time);
  DATA_VECTOR(obs_y);
  DATA_VECTOR_INDICATOR(keep,obs_y);
  DATA_STRUCT(ys_edges,directed_graph);
  DATA_STRUCT(ys_dists,dag_dists);

  DATA_IVECTOR(link_w_time);

  DATA_MATRIX(mean_design);

  DATA_IVECTOR(w_time);
  DATA_STRUCT(ws_edges,directed_graph);
  DATA_STRUCT(ws_dists,dag_dists);

  PARAMETER(mu);
  PARAMETER_VECTOR(response_pars); // Response distribution parameters except for mean
  PARAMETER_VECTOR(mean_pars); // Fixed effects B*X
  PARAMETER_VECTOR(link_w);
  PARAMETER(logtau);
  PARAMETER(logrho);
  PARAMETER(logit_w_phi);
  PARAMETER_VECTOR(proc_w);

  Type tau = exp(logtau);
  Type rho = exp(logrho);
  Type w_phi = invlogit(logit_w_phi);

  vector<vector<int> > ys_dag = ys_edges.dag;
  vector<matrix<Type> > ys_dist = ys_dists.dag_dist;

  vector<vector<int> > ws_dag = ws_edges.dag;
  vector<matrix<Type> > ws_dist = ws_dists.dag_dist;

  covariance<Type> cov(tau,rho);
  inv_link_function inv_link = {link_code};
  response_density y_density = {distribution_code};
  glm<Type> family(inv_link,
                   y_density,
                   mu,
                   mean_pars,
                   response_pars);

  vector<int> w_segment(2);
  vector<int> y_segment(2);
  vector<int> link_w_segment(2);

  vector<Type> resp_response(obs_y.size());

  Type nll = 0.0;

  // Initial time
  w_segment = get_time_segment(w_time,0);
  y_segment = get_time_segment(y_time,0);
  link_w_segment = get_time_segment(link_w_time,0);

  nngp<Type> process(cov,
                     proc_w.segment(w_segment(0),w_segment(1)),
                     0*proc_w.segment(w_segment(0),w_segment(1)), // vector of zeros
                     ws_dag,
                     ws_dist);

  observations<Type> obs(process,
                         obs_y.segment(y_segment(0),y_segment(1)),
                         keep.segment(y_segment(0),y_segment(1)),
                         ys_dag.segment(y_segment(0),y_segment(1)),
                         ys_dist.segment(y_segment(0),y_segment(1)),
                         link_w.segment(link_w_segment(0),link_w_segment(1)),
                         matrix_row_segment(mean_design,y_segment(0),y_segment(1)),
                         family);
  resp_response.segment(y_segment(0),y_segment(1)) = obs.find_response();

  nll -= process.loglikelihood()
            + obs.link_w_loglikelihood()
            + obs.y_loglikelihood();

  for(int time=1; time<n_time; time++) {
    w_segment = get_time_segment(w_time,time);
    y_segment = get_time_segment(y_time,time);
    link_w_segment = get_time_segment(link_w_time,time);

    process.update_w(proc_w.segment(w_segment(0),w_segment(1)),
                     w_phi * process.get_w());
    obs.update_y(obs_y.segment(y_segment(0),y_segment(1)),
                 keep.segment(y_segment(0),y_segment(1)),
                 ys_dag.segment(y_segment(0),y_segment(1)),
                 ys_dist.segment(y_segment(0),y_segment(1)),
                 link_w.segment(link_w_segment(0),link_w_segment(1)),
                 matrix_row_segment(mean_design,y_segment(0),y_segment(1)));
    resp_response.segment(y_segment(0),y_segment(1)) = obs.find_response();

    nll -= process.loglikelihood()
              + obs.link_w_loglikelihood()
              + obs.y_loglikelihood();
  }

  Type par_mu = mu;
  REPORT(par_mu);
  ADREPORT(par_mu);

  vector<Type> par_mean_pars = mean_pars;
  REPORT(par_mean_pars);
  ADREPORT(par_mean_pars);

  // switch statement doesn't work with REPORT and ADREPORT for some reason
  if( distribution_code == 0 ) { // Normal
    Type par_sd = exp(response_pars(0));
    REPORT(par_sd);
    ADREPORT(par_sd);
  } else if( distribution_code == 2 ) { // Neg. Binomial
    Type par_overdispersion = exp(response_pars(0))+1;
    REPORT(par_overdispersion);
    ADREPORT(par_overdispersion);
  } else {}

  Type par_tau = tau;
  REPORT(par_tau);
  ADREPORT(par_tau);
  Type working_par_logtau = logtau;
  REPORT(working_par_logtau);
  ADREPORT(working_par_logtau);

  Type par_rho = rho;
  REPORT(par_rho);
  ADREPORT(par_rho);
  Type working_par_logrho = logrho;
  REPORT(working_par_logrho);
  ADREPORT(working_par_logrho);

  Type par_w_phi = w_phi;
  REPORT(par_w_phi);
  ADREPORT(par_w_phi);
  Type working_par_logit_w_phi = logit_w_phi;
  REPORT(working_par_logit_w_phi);
  ADREPORT(working_par_logit_w_phi);

  REPORT(obs_y);

  // vector<Type> resp_response = resp_response;
  ADREPORT(resp_response);

  return(nll);
}