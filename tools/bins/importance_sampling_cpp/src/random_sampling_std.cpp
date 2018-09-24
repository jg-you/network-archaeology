/**
* \file importance_sampling.cpp
* \brief main file for importance sampling method to infere the age of edges 
in a network.
* \author Guillaume St-Onge
* \version 1.0
* \date 08/11/2017
*/

#include "DynamicNetwork.hpp"
#include "auxiliaryFunctions.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <chrono>

using namespace std;
using namespace DynNet;

typedef boost::random::uniform_01<> uniform_01;
typedef boost::random::uniform_int_distribution<> uniform_int;

const unsigned int CALIBSIZE = 100;
int main(int argc, char const *argv[])
{
	//Parameters
	string path = argv[1];
	double gamma = stod(argv[2]);
	double b = stod(argv[3]);
	unsigned int sampleSize = atoi(argv[4]);
	double seedExponent = stod(argv[5]);
	unsigned int seed;
	if (argc == 6)
	{
		seed = (unsigned int) std::chrono::high_resolution_clock::now().time_since_epoch().count();
	}
	else
	{
		seed = atoi(argv[6]);
	}

	unsigned int calibSize;
	if (sampleSize > CALIBSIZE)
	{
		calibSize = CALIBSIZE;
	}
	else
	{
		calibSize = sampleSize;
	}
	pair<double,double> param(gamma,b);

	//Load edgeList
	vector<edge> edgeList = input_edgeList(path);
	size_t M = edgeList.size();

	//Create unvarying maps associated to the observed network
	unordered_map<node,set<node> > neighborMap;
	unordered_map<node, unsigned int> degreeMap;
	edgeIntMap multiplicityMap;
	for (int i = 0; i < edgeList.size(); ++i)
	{
		//edges are assume to be in ascending order
		neighborMap[edgeList[i].first].insert(edgeList[i].second);
		neighborMap[edgeList[i].second].insert(edgeList[i].first);
		if (multiplicityMap.find(edgeList[i]) == 
			multiplicityMap.end())
		{
			multiplicityMap[edgeList[i]] = 1;			
		}
		else
		{
			multiplicityMap[edgeList[i]] += 1;
		}
		if (degreeMap.find(edgeList[i].first) != degreeMap.end())
		{
			degreeMap[edgeList[i].first] += 1;
		}
		else
		{
			degreeMap[edgeList[i].first] = 1;
		}
		if (degreeMap.find(edgeList[i].second) != degreeMap.end())
		{
			degreeMap[edgeList[i].second] += 1;
		}
		else
		{
			degreeMap[edgeList[i].second] = 1;
		}
	}
	size_t N = degreeMap.size();

	//Initialize mean marginal estimation
	unordered_map< edge, vector<double>, 
		boost::hash<edge> > meanMarginalMap;
	for (auto iter = multiplicityMap.begin(); iter != multiplicityMap.end(); 
		++iter)
	{
		vector<double> emptyVector(iter->second, 0.);
		meanMarginalMap[iter->first] = emptyVector;
	}
	double effectiveWeightSum = 0.;

	/*=======================================
		Calibration of the mean logweight
	=======================================*/

	vector<DynamicNetwork> calibSample;

	//Initialize random number generator and seed distribution
	RNGType gen(seed);
	vector<double> seedWeightVector;
	for (auto iter = edgeList.begin(); iter != edgeList.end() ; ++iter)
	{
		if (iter->first == iter->second)
		{
			seedWeightVector.push_back(0.);
		}
		else
		{
			seedWeightVector.push_back(pow(degreeMap[iter->first]
				*degreeMap[iter->second],seedExponent));
		}
	}
	boost::random::discrete_distribution<int> seedDist(seedWeightVector.begin(),
		seedWeightVector.end());

	//Get calibrating sample
	for (int n = 0; n < calibSize; ++n)
	{
		//Initialize dynamic network
		DynamicNetwork net(neighborMap, multiplicityMap);

		//choose a first random edge
		int seedIndex = seedDist(gen);
		edge seedEdge = edgeList[seedIndex];
		double seedInvWeight = 1/(seedWeightVector[seedIndex]);
		net.init(seedEdge, seedInvWeight);
		//init degree norm
		double degreeNormalization = 2.; //initial loop are impossible

		//Reconstruct the network
		while (net.get_reachableEdgeSet().size() > 0)
		{
			add_edge(net, param, degreeNormalization, gen);
		}
		calibSample.push_back(net);
	}

	//Determine mean logweight from calibration
	double meanLogweight = 0.;
	for (int n = 0; n < calibSize; ++n)
	{
		meanLogweight += calibSample[n].get_logweight();	
	}
	meanLogweight /= calibSize;

	//Estimate the mean marginal from the sample
	for (int n = 0; n < calibSize; ++n)
	{
		update_meanMarginal(calibSample[n], meanMarginalMap, meanLogweight,
			effectiveWeightSum);
	}
	calibSample.clear();

	/*=======================================
	       Refining the mean marginal
	=======================================*/

	for (int n = 0; n < sampleSize-calibSize; ++n)
	{
		//Initialize dynamic network
		DynamicNetwork net(neighborMap, multiplicityMap);

		//choose a first random edge
		int seedIndex = seedDist(gen);
		edge seedEdge = edgeList[seedIndex];
		double seedInvWeight = 1/(seedWeightVector[seedIndex]);
		net.init(seedEdge, seedInvWeight);
		//init degree norm
		double degreeNormalization = 2.; //initial loop are impossible

		//Reconstruct the network
		while (net.get_reachableEdgeSet().size() > 0)
		{
			add_edge(net, param, degreeNormalization, gen);
		}
		
		update_meanMarginal(net, meanMarginalMap, meanLogweight,
			effectiveWeightSum);
	}

	/*=======================================
	       		Output the data
	=======================================*/
	for (auto iter = meanMarginalMap.begin(); iter != meanMarginalMap.end(); 
		++iter)
	{
		double average = 0;
		for (int i = 0; i < (iter->second).size(); ++i)
		{
			average += (iter->second)[i]/effectiveWeightSum;
		}
		for (int i = 0; i < (iter->second).size(); ++i)
		{
			cout << (iter->first).first << " " << (iter->first).second << " " << i << " ";
			cout << average / double((iter->second).size()) << endl;
		}
	}

	return 0;
}