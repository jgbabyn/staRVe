#define TMB_LIB_INIT R_init_staRVe
#include <TMB.hpp>
using namespace density;

#include "include/data_in.hpp"
#include "include/time_segment.hpp"
#include "include/family.hpp"
#include "include/glm.hpp"

#include "include/covariance.hpp"
#include "include/kriging.hpp"
#include "include/nngp.hpp"
#include "include/observations.hpp"

#include "include/R_inla_barrier.hpp"

#include "model/staRVe_model.hpp"
#include "model/family.hpp"
#include "model/spde_models.hpp"



template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "staRVe_model") {
     return staRVe_model(this);
  } else if(model == "family") {
     return family(this);
  } else if(model == "barrier_model") {
    return mayo_on_sundae(this);
  } else {
     error("Unknown model.");
  }
  return Type(0.0); //Just in case..
}
