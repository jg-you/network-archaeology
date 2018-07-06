#ifndef GROWTH_MODEL_H
#define GROWTH_MODEL_H

#include <iostream>
#include <map>
#include <random>
#include "types.h"


class growth_model {
protected:
  // history of the system
  id_vec_t history_;
  id_t max_token_;
  // map of tokens to their state variables (degree, fitness, momentum, etc.)
  std::map<id_t, double_vec_t> token_states_;
  double normalization_;  // sum of all token absolute weights (weight function dependent)
  // internal distribution
  std::uniform_real_distribution<> rand_real_;
public:
  growth_model() : rand_real_(0, 1) {;}
  /// @name Virtual functions
  //{@
  // advance the model one step
  virtual void step(unsigned int t, std::mt19937& engine) {return;}
  // reset history
  virtual void intialize(unsigned int T) {history_.reserve(T); return;}
  // absolute weight of a token given its state
  virtual double weight(id_t token) const {return 1.0;}  
  //@}

  /// @name Accessors
  //@{
  id_t random_token(std::mt19937 & engine)
  {
    double accu = 0;
    double random_real = rand_real_(engine) * normalization_;
    for (auto const & kv : token_states_)
    {
      accu += weight(kv.first);
      if (random_real <= accu)
          return kv.first;
    }
    return max_token_;
  }
  void print_states(std::ostream& os) const
  {
    for (auto const & kv: token_states_)
    {
      for (auto const & v: kv.second)
      {
        os << v << " ";
      }
      os << "\n";
    }
    return;
  }
  id_vec_t get_history() const {return history_;}
  virtual void print_history(std::ostream& os) const {return;}
  //@}

  void run(unsigned int T, std::mt19937& engine)
  {
    intialize(T);
    for (unsigned int t = 0; t < T; ++t) step(t, engine);
  }
};


class unstructured_model : public growth_model {
public:
  void print_history(std::ostream& os) const
  {
    for (auto t: history_)
      os << t << "\n";
    return;
  }
};


class structured_model : public growth_model {
public:
  void print_history(std::ostream& os) const
  {
    for (unsigned int t = 0; t < history_.size(); t+=2)
    {
      if (history_[t] > history_[t + 1])
      {
        os << history_[t + 1] << " " << history_[t] << "\n";
      }
      else
      {
        os << history_[t] << " " << history_[t + 1] << "\n"; 
      }
    }
    return;
  }
};


#endif // GROWTH_MODEL_H