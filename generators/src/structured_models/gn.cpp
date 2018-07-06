#include "gn.h"

gn_model::gn_model(double gamma, unsigned int m) : structured_model()
{
  gamma_ = gamma;
  m_ = m;
}

void gn_model::intialize(unsigned int T)
{
  // clear current state
  history_.clear();
  history_.reserve(2 * T * m_);
  // single unconnected node
  token_states_[0] = {1.0};
  normalization_ = 1.0;
  max_token_ = 0;
}

void gn_model::step(unsigned int t, std::mt19937& engine)
{
  // select target nodes
  id_vec_t nodes(m_);
  for (unsigned int i = 0; i < m_; ++i)
    nodes[i] = random_token(engine);
  // add new edges and add to history
  ++max_token_;
  token_states_[max_token_] = {(double) m_};
  normalization_ += weight(max_token_); // update to normalization due to new node
  for (auto const & n: nodes)
  {
    // add to history
    history_.push_back(max_token_);
    history_.push_back(n);
    // update states
    normalization_ -= weight(n);
    ++token_states_[n][0];
    normalization_ += weight(n);
  }
}

double gn_model::weight(id_t token) const
{
  // p[0] = gamma
  return std::pow(token_states_.at(token)[0], gamma_);
}
