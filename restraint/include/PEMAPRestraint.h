/**
 *  \file IMP/isd/pEMApRestraint.h
 *  \brief A uniform/Cauchy distribution
 *  with prior knowledge on the false positive rate.
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPISD_PEMAP_RESTRAINT_H
#define IMPISD_PEMAP_RESTRAINT_H
#include "isd_config.h"
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/container_macros.h>
#include <IMP/PairContainer.h>
#include <IMP/Restraint.h>

IMPISD_BEGIN_NAMESPACE

//! Restrain atom pairs based on a set of crosslinks.
class IMPISDEXPORT pEMapRestraint : public Restraint {
public:
  //! Constructor
  /**
     \param[in] m      The IMP model 
     \param[in] c      The M*gamma normalization constant
     \param[in] psi    The xlink uncertainty nuisance particle
     \param[in] slope Exponential factor added to the score
     \param[in] part_of_log_score If true, will return probability instead of
                                  score (but the derivative will still be for a log)
     \param[in] name The restraint name

   */
  pEMapRestraint(IMP::Model* m,
                 ParticleIndex sigma,
                 Float slope = 0.01,
		 bool part_of_log_score=false,
		 std::string name = "pEMapRestraint%1%");

  //! Add a contribution to this xlink
  /**
     \param[in] ppi    The particles to which you will apply the restraint
  */

  void add_contribution(const ParticleIndex& pii0,
                        const ParticleIndex& pii1,
                        const Float& d_mic);


  //! Evaluate the restraint just for a subset of contribution indexes
  Float evaluate_for_contributions(Ints c, DerivativeAccumulator *accum) const;

  //! Get number of contributions added to the restraint
  unsigned int get_number_of_contributions() const { return ppis_.size(); }

  //! Get the particles in a contribution
  ParticleIndexPair get_contribution(int i) const { return ppis_[i]; }

  //! Get the psi nuisance in a contribution
  //ParticleIndex get_sigma() const { return sigma_; }

  virtual double unprotected_evaluate(IMP::DerivativeAccumulator* accum) const IMP_OVERRIDE;
  virtual IMP::ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE;
  void show(std::ostream &out) const { out << "pEMap restraint with "
                                           <<get_number_of_contributions()
                                           <<" contributions"; }
  void set_slope(Float slope) { slope_=slope; }
  void set_part_of_log_score(bool hey) { part_of_log_score_=hey; }

  IMP_OBJECT_METHODS(pEMapRestraint);
 private:
  ParticleIndexPairs ppis_;
  ParticleIndex sigma_;
  Float slope_;
  bool part_of_log_score_;
  std::vector<double> d_mics_;
  Ints default_range_;
};

IMPISD_END_NAMESPACE

#endif /* IMPISD_PEMAP_RESTRAINT_H */
