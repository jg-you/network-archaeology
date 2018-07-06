#include "gn.h"

generalized_gn_model::generalized_gn_model(double a,
                               double alpha,
                               double b,
                               unsigned int m1,
                               unsigned int m2,
                               double gamma,
                               double tau,
                               bool directed) : structured_model()
{
  a_ = a;     
  alpha_ = alpha;     
  b_ = b;     
  m1_ = m1;     
  m2_ = m2;     
  gamma_ = gamma;     
  tau_ = tau;     
  directed_ = directed;
}

void generalized_gn_model::intialize(unsigned int T)
{
  // clear current state
  history_.clear();
  history_.reserve(2 * T * (m1_ + m2_));
  adj_list_.clear();
  adj_list_.reserve(T * (m1_ + m2_));
  // adj_list_.push_back(neighbourhood_t());
  // unconnected node
  token_states_[0] = {1.0};
  token_states_[1] = {1.0};
  normalization_ = 2.0;
  max_token_ = 1;
  adj_list_.push_back({1});
  adj_list_.push_back({0});
  history_.push_back(0);
  history_.push_back(1);
}

void generalized_gn_model::step(unsigned int t, std::mt19937& engine)
{
  if (rand_real_(engine) < p(t)) // new node
  {
    // select target nodes
    id_vec_t nodes(m1_);
    for (unsigned int i = 0; i < m1_; ++i)
      nodes[i] = random_token(engine);
    // add new edges and add to history
    for (auto const & n: nodes)
    {
      ++max_token_;
      // add to history
      history_.push_back(max_token_);
      history_.push_back(n);
      // update states
      token_states_[max_token_] = {1.0};
      normalization_ += 1.0; // update to normalization due to new node
      normalization_ -= weight(n);
      ++token_states_[n][0];
      normalization_ += weight(n);
      // connect
      adj_list_.push_back({n});
      if (!directed_)
        adj_list_[n].insert(max_token_);
    }
  }
  else  // densification
  {
    id_vec_t source(m2_);
    for (unsigned int i = 0; i < m2_; ++i)
    {
      source[i] = random_token(engine);
    }

    for (unsigned int i = 0; i < m2_; ++i)
    {
      id_t target;
      target = random_token(engine);
      history_.push_back(source[i]);
      history_.push_back(target);
      if (source[i] != target)
      {
          normalization_ -= weight(source[i]);
          normalization_ -= weight(target);
          ++token_states_[source[i]][0];
          ++token_states_[target][0];
          normalization_ += weight(source[i]);
          normalization_ += weight(target);
      }
      else  // same node
      {
          normalization_ -= weight(target);
          token_states_[target][0] += 2;
          normalization_ += weight(target);        
      }
      adj_list_[source[i]].insert(target);
      if (!directed_)
        adj_list_[target].insert(source[i]);
    }
  }
}

double generalized_gn_model::weight(id_t token) const
{
  // p[0] = gamma
  return std::pow(token_states_.at(token)[0], gamma_);
}

double generalized_gn_model::p(unsigned int t) const
{
  return a_ * std::pow((double) t + tau_, - alpha_) + b_;
}
