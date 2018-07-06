#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <set>
#include <utility>

typedef unsigned int id_t;

typedef std::pair<id_t, id_t> edge_t;
typedef std::vector<edge_t> edge_list_t;
typedef std::multiset<id_t> neighbourhood_t;
typedef std::vector<neighbourhood_t> adj_list_t;

typedef struct mcmc_move_t
{
  id_t t;
  id_t token;
} mcmc_move_t;

typedef std::vector<unsigned int> uint_vec_t;
typedef std::vector<id_t> id_vec_t;
typedef std::vector<int> int_vec_t;
typedef std::vector<float> float_vec_t;
typedef std::vector<double> double_vec_t;
typedef std::vector< std::vector<unsigned int> > uint_mat_t;
typedef std::vector< std::vector<int> > int_mat_t;
typedef std::vector< std::vector<float> > float_mat_t;

#endif // TYPES_H
