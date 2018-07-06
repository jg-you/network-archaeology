#ifndef GN_H
#define GN_H

#include <cmath>

#include "../growth_model.h"
#include "../types.h"

class gn_model : public structured_model {
public:
  // no seed constructor
  gn_model(double gamma, unsigned int m);
  void intialize(unsigned int T);
  // overloaded step function
  void step(unsigned int t, std::mt19937& engine);
  // model specific
  double weight(id_t token) const;

private:
   double gamma_;
   unsigned int m_;
};

class generalized_gn_model : public structured_model {
public:
  // no seed constructor
  generalized_gn_model(double a, double alpha, double b, unsigned int m1, unsigned int m2, double gamma, double tau, bool directed=false);
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
  double gamma_;
  double tau_;
  adj_list_t adj_list_;
  bool directed_;
};


#endif // GN_H