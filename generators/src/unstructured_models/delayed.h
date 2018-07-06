#ifndef DELAYED_H
#define DELAYED_H

#include <cmath>

#include "../growth_model.h"
#include "../types.h"

class delayed_model : public unstructured_model {
public:
  // no seed constructor
  delayed_model(double a, double alpha, double b, double mu, double tau);
  void intialize(unsigned int T);
  // overloaded step function
  void step(unsigned int t, std::mt19937& engine);
  // model specific
  double p(unsigned int t) const;
  double weight(id_t token) const;

private:
  double_vec_t params_;
};

#endif // DELAYED_H