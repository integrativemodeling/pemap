/**
 *  \file isd/pEMapRestraint.h
 *  \brief A sigmoid shaped restraint between
 *  residues with discrete classifier
 *  and ambiguous assignment. To be used with
 *  cross-linking mass-spectrometry data.
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/isd/pEMapRestraint.h>
#include <IMP/algebra/VectorD.h>
#include <IMP/core/XYZ.h>
#include <IMP/isd/Scale.h>
#include <boost/math/special_functions/erf.hpp>
#include <math.h>

IMPISD_BEGIN_NAMESPACE

namespace {
  static const double PI = 3.1415926535897931;
}
pEMapRestraint::pEMapRestraint(Model* m,
                               ParticleIndex sigma,
                               Float slope,
			       bool part_of_log_score,
			       std::string name):
  
  Restraint(m, name),
  sigma_(sigma),
  slope_(slope),
  part_of_log_score_(part_of_log_score)
  {
}

void pEMapRestraint::add_contribution(const ParticleIndex& pii0,
                                      const ParticleIndex& pii1,
				      const Float& d_mic){
  
  ppis_.push_back(ParticleIndexPair(pii0,pii1));
  d_mics_.push_back(d_mic);
  default_range_.push_back((int)default_range_.size());
}

Float pEMapRestraint::evaluate_for_contributions(Ints c,
                                                 DerivativeAccumulator *accum) const{

  double score = 0.0;
  
  Float sigma = isd::Scale(get_model(),sigma_).get_scale();
  
  // Loop over the contributions and do scoring 
  for (Ints::const_iterator nit=c.begin();nit!=c.end();++nit){
    int n = *nit;
    core::XYZ d0(get_model(),ppis_[n][0]);
    core::XYZ d1(get_model(),ppis_[n][1]);

    Float dm = algebra::get_distance(d0.get_coordinates(),
                                     d1.get_coordinates());
    
    if (dm > d_mics_[n]){
      double prob = std::exp(-(dm-d_mics_[n])*(dm-d_mics_[n])/(2*sigma*sigma));
      //double prior = std::exp(-(dm-d_mics_[n])*(dm-d_mics_[n])/(2*sigma*sigma));
      score += -std::log(prob);
      //std::cout<<"h: "<<dm<<" "<<d_mics_[n]<<" "<<-std::log(prob)<<std::endl;
    }
  }
  return score;
}

double pEMapRestraint::unprotected_evaluate(DerivativeAccumulator *accum)
  const {
  return evaluate_for_contributions(default_range_,accum);
}

ModelObjectsTemp pEMapRestraint::do_get_inputs() const {
    ParticlesTemp ret;
    for (unsigned int k = 0; k < get_number_of_contributions(); ++k) {
        ret.push_back(get_model()->get_particle(ppis_[k][0]));
        ret.push_back(get_model()->get_particle(ppis_[k][1]));
    }
    return ret;
}

IMPISD_END_NAMESPACE
