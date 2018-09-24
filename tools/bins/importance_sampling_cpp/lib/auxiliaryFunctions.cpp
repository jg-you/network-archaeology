/**
* \file auxiliaryFunctions.cpp
* \brief functions for the module auxiliaryFunctions
* \author Guillaume St-Onge
* \version 1.0
* \date 08/11/2017
*/

#include "auxiliaryFunctions.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>

using namespace std;

namespace DynNet
{//start of namespace DynNet

typedef boost::random::uniform_int_distribution<> uniform_int;


/**
 * \fn void order(edge p_edge)
 * \brief Assure the edge is in ascending order
 * \param[in] p_edge 
 */
void order(edge& p_edge)
{
	if (p_edge.first > p_edge.second)
	{
		swap(p_edge.first,p_edge.second);
	}
}

/**
 * \fn vector<edge> input_edgeList(string p_path)
 * \brief Input the edge list from file
 * \param[in] p_path path name to the file
 */
vector<edge> input_edgeList(string p_path)
{
	ifstream inStream;
	inStream.open(p_path, ios::in);
	string line;

	vector<edge> edgeList;

	while(getline(inStream,line))
	{
		edge e;
		stringstream line_stream(line);
		line_stream >> e.first >> e.second;
		//swap edge to assure ascending order
		order(e);
		edgeList.push_back(e);
	}
	inStream.close();

	return edgeList;
}

/**
 * \fn pair<double,double> update_estimator(DynamicNetwork& p_net, edge& p_chosenEdge,
	pair<double,double>& p_est_param)ble,double>& p_model_param,
	edge& p_chosenEdge, double& p_degreeNormalization)
 * \brief Get likelihood for the new edge and update normalisation.
 * \param[in] p_net Structure of a network in reconstruction.
 * \param[in] p_chosenEdge new edge
 * \param[in] p_est_param Estimated parameter
 */
pair<double,double> update_estimator(DynamicNetwork& p_net, edge& p_chosenEdge,
	pair<double,double>& p_est_param)
{
	double n = p_net.get_currentDegreeMap().size();
	double m = p_net.get_currentEdgeList().size();
	double kavg = 2*m/n;
	double ksum = kavg*n;
	double k2avg = p_est_param.first + kavg*kavg;
	double k2sum = n*k2avg;

	unordered_map<node, unsigned int> currentDegreeMap = p_net.get_currentDegreeMap();

	//update parameter according to new node
	m += 1;
	ksum += 2;
	if (currentDegreeMap.find(p_chosenEdge.first) != 
		currentDegreeMap.end())
	{
		if (currentDegreeMap.find(p_chosenEdge.second) != 
		currentDegreeMap.end())
		{
			//2 old nodes, n stays the same
			//check if self-loop
			if (p_chosenEdge.first == p_chosenEdge.second)
			{
				k2sum -= (currentDegreeMap[p_chosenEdge.first]) * (currentDegreeMap[p_chosenEdge.first]);
				k2sum += (currentDegreeMap[p_chosenEdge.first] + 2) * (currentDegreeMap[p_chosenEdge.first] + 2);
				// k2sum -= pow(currentDegreeMap[p_chosenEdge.first],2.);
				// k2sum += pow(currentDegreeMap[p_chosenEdge.first] + 2,2.);
			}
			else
			{
			// 	k2sum -= pow(currentDegreeMap[p_chosenEdge.first],2.);
			// 	k2sum += pow(currentDegreeMap[p_chosenEdge.first] + 1,2.);
			// 	k2sum -= pow(currentDegreeMap[p_chosenEdge.second],2.);
			// 	k2sum += pow(currentDegreeMap[p_chosenEdge.second] + 1,2.);
				k2sum -= (currentDegreeMap[p_chosenEdge.first]) * (currentDegreeMap[p_chosenEdge.first]);
				k2sum += (currentDegreeMap[p_chosenEdge.first] + 1) * (currentDegreeMap[p_chosenEdge.first] + 1);
				k2sum -= (currentDegreeMap[p_chosenEdge.second]) * (currentDegreeMap[p_chosenEdge.second]);
				k2sum += (currentDegreeMap[p_chosenEdge.second] + 1) * (currentDegreeMap[p_chosenEdge.second] + 1);
			}
		}
		else
		{
			//first is old, second is new
			n += 1;
			// k2sum -= pow(currentDegreeMap[p_chosenEdge.first],2.);
			// k2sum += pow(currentDegreeMap[p_chosenEdge.first] + 1,2.);
			k2sum -= (currentDegreeMap[p_chosenEdge.first]) *  (currentDegreeMap[p_chosenEdge.first]);
			k2sum += (currentDegreeMap[p_chosenEdge.first] + 1) * (currentDegreeMap[p_chosenEdge.first] + 1);
			k2sum += 1;
		}
	}
	else
	{
		//first is new, second is old
		n += 1;
		// k2sum -= pow(currentDegreeMap[p_chosenEdge.second],2.);
		// k2sum += pow(currentDegreeMap[p_chosenEdge.second] + 1,2.);
		k2sum -= (currentDegreeMap[p_chosenEdge.second]) * (currentDegreeMap[p_chosenEdge.second]);
		k2sum += (currentDegreeMap[p_chosenEdge.second] + 1) * (currentDegreeMap[p_chosenEdge.second] + 1);
		k2sum += 1;
	}

	// pair<double,double> new_est_param(k2sum/n - pow(ksum/n, 2.),(n-2)/(m-1));
	pair<double,double> new_est_param(k2sum/n - (ksum/n) * (ksum/n), (n-2)/(m-1));

	return new_est_param;
}
		

/**
 * \fn double update_modelLikelihood(DynamicNetwork& p_net, pair<double,
 	double>& p_model_param,	edge& p_chosenEdge, double& p_degreeNormalization)
 * \brief Get likelihood for the new edge and update normalisation.
 * \param[in] p_net Structure of a network in reconstruction.
 * \param[in] p_model_param Model parameters
 * \param[in] p_chosenEdge new edge
 * \param[in] p_degreeNormalization Normalization factor for probability
 */
double update_modelLikelihood(DynamicNetwork& p_net, 
	pair<double,double>& p_model_param, edge& p_chosenEdge, 
	double& p_degreeNormalization)
{
	double model_likelihood = 1.;
	if (p_net.get_currentDegreeMap().find(p_chosenEdge.first) != 
		p_net.get_currentDegreeMap().end())
	{
		if (p_net.get_currentDegreeMap().find(p_chosenEdge.second) != 
		p_net.get_currentDegreeMap().end())
		{
			//2 old nodes
			double k1 = p_net.get_currentDegreeMap().at(p_chosenEdge.first);
			double k2 = p_net.get_currentDegreeMap().at(p_chosenEdge.second);
			model_likelihood *= (1 - p_model_param.second)*pow(k1*k2,
				p_model_param.first)/pow(p_degreeNormalization,2.);
			//check if self-loop
			if (p_chosenEdge.first == p_chosenEdge.second)
			{
				p_degreeNormalization -= pow(k1,p_model_param.first);
				p_degreeNormalization += pow(k1+2,p_model_param.first);
			}
			else
			{
				p_degreeNormalization -= pow(k1,p_model_param.first);
				p_degreeNormalization -= pow(k2,p_model_param.first);
				p_degreeNormalization += pow(k1+1,p_model_param.first);
				p_degreeNormalization += pow(k2+1,p_model_param.first);
			}
		}
		else
		{
			//first is old, second is new
			double k1 = p_net.get_currentDegreeMap().at(p_chosenEdge.first);
			model_likelihood *= p_model_param.second*pow(k1,
				p_model_param.first)/p_degreeNormalization;
			p_degreeNormalization -= pow(k1,p_model_param.first);
			p_degreeNormalization += pow(k1+1,p_model_param.first)+1;
		}
	}
	else
	{
		//first is new, second is old
		double k2 = p_net.get_currentDegreeMap().at(p_chosenEdge.second);
		model_likelihood *= p_model_param.second*pow(k2,
			p_model_param.first)/p_degreeNormalization;
		p_degreeNormalization -= pow(k2,p_model_param.first);
		p_degreeNormalization += pow(k2+1,p_model_param.first)+1;
	}

	return model_likelihood;
}

/**
 * \fn void add_edge(DynamicNetwork& p_net, pair<double,double>& p_model_param,
	pair<double,double>& p_target_param, pair<double,double>& p_est_param, 
	pair<double,double>& p_bias_param, RNGType& gen)
 * \brief Add a new edge to the network following to the growth model.
 * \param[in] p_net Structure of a network in reconstruction.
 * \param[in] p_model_param Model parameters
 * \param[in] p_target_param Target parameter for bias
 * \param[in] p_est_param Estimated parameter for bias
 * \param[in] p_bias_param Tunable parameter for bias
 * \param[in] p_degreeNormalization Normalization factor for probability
 * \param[in] gen Random number generator.
 */
void add_edge(DynamicNetwork& p_net, pair<double,double>& p_model_param,
	pair<double,double>& p_target_param, pair<double,double>& p_est_param, 
	pair<double,double>& p_bias_param, double& p_degreeNormalization,
	RNGType& gen)
{
	//this add_edge method is used with important_sampling

	edgeSet reachableEdgeSet = p_net.get_reachableEdgeSet();
	unordered_map<unsigned int, edge> indexMap;
	vector<double> weightList(reachableEdgeSet.size(), 0);
	// vector<double> weightList;
	double weightNorm = 0.;
	bool oldSeed = false;
	bool oldTarget = false;

	//get weights
	unsigned int index = 0;
	for (auto iter = reachableEdgeSet.begin(); 
		iter != reachableEdgeSet.end(); ++iter)
	{
		indexMap[index] = *iter;

		double weight = p_net.get_remainingEdgeMap().at(*iter); //bias for multiple edges
		edge possibleEdge = *iter;
		pair<double,double> new_est_param = update_estimator(p_net, possibleEdge, 
			p_est_param);
		// weight *= exp(-p_bias_param.first*(new_est_param.second 
		// 	- p_target_param.second));
		// weight *= exp(-p_bias_param.second*(new_est_param.first
		// 	- p_target_param.first));
		weight *= exp(-p_bias_param.first*(new_est_param.second- p_target_param.second) * (new_est_param.second- p_target_param.second));  // density
		weight *= exp(-p_bias_param.second*(new_est_param.first- p_target_param.first) * (new_est_param.first- p_target_param.first));  // variance

		// weightList.push_back(weight);

		weightList[index] = weight;
		weightNorm += weight;

		index += 1;
	}

	//Get a new random edge
	boost::random::discrete_distribution<int> edgeDist(weightList.begin(), weightList.end());
	unsigned int newEdgeIndex = edgeDist(gen);
	edge chosenEdge = indexMap[newEdgeIndex];

	//determine the model likelihood
	double model_likelihood = update_modelLikelihood(p_net, p_model_param, chosenEdge, 
		p_degreeNormalization);

	//update the estimator
	p_est_param = update_estimator(p_net, chosenEdge, p_est_param);
	
	// cout << p_est_param.first << " " << p_est_param.second << "\n";

	//add the edge with according probability ratio (model_likelihood/bias_prob)
	double ratio = model_likelihood*weightNorm/(weightList[newEdgeIndex]);
	//account for multiple edge indistinguishability
	ratio *= p_net.get_remainingEdgeMap().at(chosenEdge);

	p_net.add(chosenEdge, ratio);
}

/**
 * \fn void add_edge(DynamicNetwork& p_net, pair<double,double>& p_param,
	double& p_degreeNormalization, RNGType& gen)
 * \brief Add a new edge to the network following to the growth model.
 * \param[in] p_net Structure of a network in reconstruction.
 * \param[in] p_param Model parameters
 * \param[in] p_degreeNormalization Normalization factor for probability
 * \param[in] gen Random number generator.
 */
void add_edge(DynamicNetwork& p_net, pair<double,double>& p_param,
	double& p_degreeNormalization, RNGType& gen)
{
	edgeSet reachableEdgeSet = p_net.get_reachableEdgeSet();
	uniform_int indexDist(0, p_net.get_reachableEdge()-1);

	// cout << reachableEdgeSet.size() << " " << p_net.get_reachableEdge() << endl;
	int index = indexDist(gen);

	auto iter = reachableEdgeSet.begin();
	//count multiple edge with remaining map
	while (index >= 0)
	{
		unsigned int remaining = p_net.get_remainingEdgeMap().at(*iter);
		index -= remaining;
		if (index >= 0)
		{
			iter++ ;
		}
	}

	bool oldSeed = false;
	bool oldTarget = false;
	if (p_net.get_currentDegreeMap().find(iter->first) != 
		p_net.get_currentDegreeMap().end())
	{
		oldSeed = true;
	}
	else
	{
		oldSeed = false;
	}
	if (p_net.get_currentDegreeMap().find(iter->second) != 
		p_net.get_currentDegreeMap().end())
	{
		oldTarget = true;
	}
	else
	{
		oldTarget = false;
	}

	double weight = p_net.get_reachableEdge();
	if (oldSeed)
	{
		if (oldTarget)
		{
			double k1 = p_net.get_currentDegreeMap().at(iter->first);
			double k2 = p_net.get_currentDegreeMap().at(iter->second);
			weight *= (1- p_param.second)/pow(p_degreeNormalization,2.);
			weight *= pow(k1, p_param.first);
			weight *= pow(k2, p_param.first);
			if (iter->second == iter->first)	
			{
				p_degreeNormalization -= pow(k1, p_param.first);
				p_degreeNormalization += pow(k1+2, p_param.first);
			}
			else
			{
				p_degreeNormalization -= pow(k1, p_param.first);
				p_degreeNormalization -= pow(k2, p_param.first);
				p_degreeNormalization += pow(k1+1., p_param.first);
				p_degreeNormalization += pow(k2+1., p_param.first);
			}
		}
		else
		{
			double k1 = p_net.get_currentDegreeMap().at(iter->first);
			weight *= p_param.second/p_degreeNormalization;
			weight *= pow(k1, p_param.first);
			p_degreeNormalization -= pow(k1, p_param.first);
			p_degreeNormalization += pow(k1+1., p_param.first) + 1;
		}
	}
	else
	{
		if (oldTarget)
		{
			double k2 = p_net.get_currentDegreeMap().at(iter->second);
			weight *= p_param.second/p_degreeNormalization;
			weight *= pow(k2, p_param.first);
			// weight *= p_net.get_remainingEdgeMap().at(*iter);//bias to account
			p_degreeNormalization -= pow(k2, p_param.first);
			p_degreeNormalization += pow(k2+1., p_param.first) + 1;
		}
		else
		{
			//impossible to have 2 new nodes
			cout << "ERROR : two new nodes" << endl;
		}
	}

	p_net.add(*iter, weight);
}

/**
 * \fn void update_meanMarginal(DynamicNetwork& p_net, unordered_map< edge, 
 	vector<double>, boost::hash<edge> > p_meanMarginalMap, 
 	double meanLogweight, double& p_effectiveWeightSum)
 * \brief Update mean marginal evaluation from a network sample.
 * \param[in] p_net Structure of a network in reconstruction.
 * \param[in] p_meanMarginalMap Map of the mean marginal.
 * \param[in] p_meanLogweight Mean value of the logweight.
 * \param[in] p_effectiveWeightSum Normalization value to update.
 */
void update_meanMarginal(DynamicNetwork& p_net, unordered_map< edge, 
	vector<double>, boost::hash<edge> >& p_meanMarginalMap, 
	double meanLogweight, double& p_effectiveWeightSum)
{
	edgeIntMap edgeCountMap;
	double effectiveWeight = exp(p_net.get_logweight()-meanLogweight);
	double t = 0.;
	for (auto iter = p_net.get_currentEdgeList().begin(); 
		iter != p_net.get_currentEdgeList().end(); ++iter)
	{
		p_meanMarginalMap[*iter][edgeCountMap[*iter]] += t*effectiveWeight;
		edgeCountMap[*iter] += 1;
		t += 1.;
	}
	p_effectiveWeightSum += effectiveWeight;
}

/**
 * \fn bool update_meanMarginal(DynamicNetwork& p_net, unordered_map< edge, 
 	vector<double>, boost::hash<edge> > p_meanMarginalMap,
 	double meanLogweight, double& p_effectiveWeightSum, double p_rtol)
 * \brief Update mean marginal evaluation from a network sample.
 * \param[in] p_net Structure of a network in reconstruction.
 * \param[in] p_meanMarginalMap Map of the mean marginal.
 * \param[in] p_meanLogweight Mean value of the logweight.
 * \param[in] p_effectiveWeightSum Normalization value to update.
 * \param[in] p_rtol Tolerance for convergence.
 */
bool update_meanMarginal(DynamicNetwork& p_net, unordered_map< edge, 
	vector<double>, boost::hash<edge> >& p_meanMarginalMap, 
	double meanLogweight, double& p_effectiveWeightSum, double p_rtol,
	unsigned int sampleSize)
{
	edgeIntMap edgeCountMap;
	bool converged = true;
	double effectiveWeight = exp(p_net.get_logweight()-meanLogweight);
	p_effectiveWeightSum += effectiveWeight;
	double t = 0.;
	for (auto iter = p_net.get_currentEdgeList().begin(); 
		iter != p_net.get_currentEdgeList().end(); ++iter)
	{
		p_meanMarginalMap[*iter][edgeCountMap[*iter]] += t*effectiveWeight;
		double avg_t = (p_meanMarginalMap[*iter][edgeCountMap[*iter]]
			/p_effectiveWeightSum);
		if (abs(t-avg_t)/sqrt(sampleSize) > p_rtol) //comportement moyen
		{
			converged = false;
		}
		edgeCountMap[*iter] += 1;
		t += 1.;
	}

	return converged;
}

}//end of namespace DynNet