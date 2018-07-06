#include "delayed.h"

delayed_model::delayed_model(double a, double alpha, double b, double mu, double tau) : unstructured_model()
{
  params_ = {a, alpha, b, mu, tau};
}

void delayed_model::intialize(unsigned int T)
{
  // clear current state
  token_states_.clear();
  history_.clear();
  history_.reserve(T);
  // add a single token
  history_.push_back(0);
  max_token_ = 0;
  token_states_[0] = {1.0};
  normalization_ = 1.0;
}

void delayed_model::step(unsigned int t, std::mt19937& engine)
{
  if (rand_real_(engine) < p(t))
  {
    ++max_token_;
    history_.push_back(max_token_);
    token_states_[max_token_] = {1.0};
    normalization_ += 1.0;
  }
  else
  { 
    id_t token = random_token(engine);
    history_.push_back(token);
    normalization_ -= weight(token); // remove old weight from normalization.
    ++token_states_[token][0]; // token now appears one more time.
    normalization_ += weight(token); // add new weight.
  }
}

double delayed_model::p(unsigned int t) const
{
  return params_[0] * std::pow((double) t + params_[4], -params_[1]) + params_[2];
}

double delayed_model::weight(id_t token) const
{
  return std::pow(token_states_.at(token)[0], params_[3]);
}
