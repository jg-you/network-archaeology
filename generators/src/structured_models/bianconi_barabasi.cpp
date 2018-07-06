#include "bianconi_barabasi.h"

bianconi_barabasi_model::bianconi_barabasi_model(double mu, double nu, double avg, double std, unsigned int m)
: structured_model(), normal_distribution_(avg, std)
{
  mu_ = mu;
  nu_ = nu;
  m_ = m;
  avg_ = avg;
}

void bianconi_barabasi_model::intialize(unsigned int T)
{
  // clear current state
  history_.clear();
  history_.reserve(2 * T * m_);
  // single unconnected node
  token_states_[0] = {1.0, avg_};
  normalization_ = weight(0);
  max_token_ = 0;
}

void bianconi_barabasi_model::step(unsigned int t, std::mt19937& engine)
{
  // select target nodes
  id_vec_t nodes(m_);
  for (unsigned int i = 0; i < m_; ++i)
    nodes[i] = random_token(engine);
  // add new edges and add to history
  for (auto const & n: nodes)
  {
    ++max_token_;
    // add to history
    history_.push_back(max_token_);
    history_.push_back(n);
    // update states
    token_states_[max_token_] = {1.0, normal_distribution_(engine)};
    normalization_ += weight(max_token_); // update to normalization due to new node
    normalization_ -= weight(n);
    ++token_states_[n][0];
    normalization_ += weight(n);
  }
}

double bianconi_barabasi_model::weight(id_t token) const
{
  // p[0] = mu
  auto state = token_states_.at(token);
  return std::pow(state[0], mu_ ) * std::pow(state[1], nu_);
}
