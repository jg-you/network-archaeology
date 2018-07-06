#ifndef GN_GEN_MODEL_H
#define GN_GEN_MODEL_H

#include <cmath>

#include "growth_model.h"
#include "types.h"

class gn_gen_model_t : public structured_model_t {
public:
  // no seed constructor
  gn_gen_model_t(double a, double alpha, double b, unsigned int m1, unsigned int m2, double mu, double tau, bool directed=false);
  void intialize(unsigned int T);
  // overloaded step function
  void step(unsigned int t, std::mt19937& engine);
  // model specific
  double weight(id_t token) const;
  double p(unsigned int t) const;

private:
  double a_;
  double alpha_;
  double b_;
  unsigned int m1_;
  unsigned int m2_;
  double mu_;
  double tau_;
  adj_list_t adj_list_;
  bool directed_;
};

#endif // GN_GEN_MODEL_H