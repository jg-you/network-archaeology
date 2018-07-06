#ifndef BIANCONI_BARABASI_H
#define BIANCONI_BARABASI_H

#include <cmath>
#include <random>

#include "../growth_model.h"
#include "../types.h"

class bianconi_barabasi_model : public structured_model {
public:
  // no seed constructor
  bianconi_barabasi_model(double mu, double nu, double avg=1, double std=0, unsigned int m=1);
  void intialize(unsigned int T);
  // overloaded step function
  void step(unsigned int t, std::mt19937& engine);
  // model specific
  double weight(id_t token) const;

private:
  double mu_;
  double nu_;
  unsigned int m_;
  double avg_;
  std::normal_distribution<double> normal_distribution_;
};

#endif // BIANCONI_BARABASI_H