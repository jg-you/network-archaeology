/**
* \file auxiliaryFunctions.hpp
* \brief Module for auxiliary functions
* \author Guillaume St-Onge
* \version 1.0
* \date 08/11/2017
*/

#include "DynamicNetwork.hpp"
#include <string>
#include <boost/random/mersenne_twister.hpp>

#ifndef AUXILIARYFUNCTIONS_HPP_
#define AUXILIARYFUNCTIONS_HPP_

namespace DynNet
{//start of namespace DynNet

typedef boost::mt19937 RNGType;


void order(edge& p_edge);
std::vector<edge> input_edgeList(std::string p_path);
std::pair<double,double> update_estimator(DynamicNetwork& p_net, edge& p_chosenEdge,
	std::pair<double,double>& p_est_param);
double update_modelLikelihood(DynamicNetwork& p_net, 
	std::pair<double,double>& p_model_param,
	edge& p_chosenEdge, double& p_degreeNormalization);
void add_edge(DynamicNetwork& p_net, std::pair<double,double>& p_model_param,
	std::pair<double,double>& p_target_param, 
	std::pair<double,double>& p_est_param, 
	std::pair<double,double>& p_bias_param, 
	double& p_degreeNormalization, RNGType& gen);
void add_edge(DynamicNetwork& p_net, std::pair<double,double>& p_param,
	double& p_degreeNormalization, RNGType& gen);

void update_meanMarginal(DynamicNetwork& p_net, std::unordered_map< edge, 
	std::vector<double>, boost::hash<edge> >& p_meanMarginalMap, 
	double meanLogweight, double& p_effectiveWeightSum);

bool update_meanMarginal(DynamicNetwork& p_net, std::unordered_map< edge, 
	std::vector<double>, boost::hash<edge> >& p_meanMarginalMap, 
	double meanLogweight, double& p_effectiveWeightSum, double p_rtol,
	unsigned int sampleSize);

}//end of namespace DynNet

#endif /* AUXILIARYFUNCTIONS_HPP_ */