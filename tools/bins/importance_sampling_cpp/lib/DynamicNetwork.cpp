/**
* \file DynamicNetwork.cpp
* \brief Methods for the class DynamicNetwork
* \author Guillaume St-Onge
* \version 1.0
* \date 08/11/2017
*/

#include "DynamicNetwork.hpp"
#include "auxiliaryFunctions.hpp"
#include <cmath>

using namespace std;

namespace DynNet
{//start of namespace DynNet


/*---------------------------
 *      Constructor
 *---------------------------*/

/**
* \brief Constructor of the class
* \param[in] 
*/
DynamicNetwork::DynamicNetwork(unordered_map<node,set<node> >& p_neighborMap,
    	edgeIntMap& p_multiplicityMap) : 
	m_neighborMap(p_neighborMap), m_currentEdgeList(), m_currentEdgeSet(),
	m_reachableEdgeSet(), m_currentDegreeMap(), 
	m_remainingEdgeMap(p_multiplicityMap), m_logweight(0.), m_reachableEdge(0)
{
}


/*---------------------------
 *         Mutators
 *---------------------------*/

/**
* \brief Initialize the network with a seed
* \param[in] p_seedEdge Edge which is the seed of the growing network
* \param[in] p_seedInvWeight Factor to account for the marginal
*/
void DynamicNetwork::init(edge p_seedEdge, double p_seedInvWeight)
{
	m_currentEdgeSet.insert(p_seedEdge);
	m_currentEdgeList.push_back(p_seedEdge);
	m_logweight += log(p_seedInvWeight);
	update_reachableEdge(p_seedEdge);
}

/**
* \brief Initialize the network with a seed
* \param[in] p_seedEdge Edge which is the seed of the growing network
*/
void DynamicNetwork::init(edge p_seedEdge)
{
	m_currentEdgeSet.insert(p_seedEdge);
	m_currentEdgeList.push_back(p_seedEdge);
	update_reachableEdge(p_seedEdge);
}

/**
* \brief Add an edge to the network
* \param[in] p_edge Edge to add to the network
* \param[in] p_weightNorm Normalization factor for the edges at that time
*/
void DynamicNetwork::add(edge p_edge, double p_weightNorm)
{
	m_currentEdgeSet.insert(p_edge);
	m_currentEdgeList.push_back(p_edge);
	m_logweight += log(p_weightNorm);
	update_reachableEdge(p_edge); //newly reachable edges
}

/**
 * \brief Update the set of edges reachable for reconstruction by the edge
 * \param[in] p_edge Edge used to determine reachable edges
 */
void DynamicNetwork::update_reachableEdge(edge p_edge)
{
	//get new reachable edges
	for (auto iter = m_neighborMap.at(p_edge.first).begin(); 
		iter != m_neighborMap.at(p_edge.first).end(); ++iter)
	{
		edge e(p_edge.first,*iter);
		order(e);
		if (m_reachableEdgeSet.count(e) == 0 and m_remainingEdgeMap.at(e) > 0)
		{	//new edges available
			m_reachableEdge += m_remainingEdgeMap.at(e);
			m_reachableEdgeSet.insert(e);
		}
	}
	for (auto iter = m_neighborMap.at(p_edge.second).begin(); 
		iter != m_neighborMap.at(p_edge.second).end(); ++iter)
	{
		edge e(p_edge.second,*iter);
		order(e);
		if (m_reachableEdgeSet.count(e) == 0 and m_remainingEdgeMap.at(e) > 0)
		{	//new edges available
			m_reachableEdge += m_remainingEdgeMap.at(e);
			m_reachableEdgeSet.insert(e);
		}
	}

	//remove the new edge
	m_remainingEdgeMap[p_edge] -= 1;
	m_reachableEdge -= 1;
	if (m_remainingEdgeMap[p_edge] == 0)
	{
		m_reachableEdgeSet.erase(p_edge);
	}

	//update the degree of nodes
	if (m_currentDegreeMap.find(p_edge.first) != m_currentDegreeMap.end())
	{
		m_currentDegreeMap[p_edge.first] += 1;
	}
	else
	{
		m_currentDegreeMap[p_edge.first] = 1;
	}
	if (m_currentDegreeMap.find(p_edge.second) != m_currentDegreeMap.end())
	{
		m_currentDegreeMap[p_edge.second] += 1;
	}
	else
	{
		m_currentDegreeMap[p_edge.second] = 1;
	}
}

/*---------------------------
 *         Accessor
 *---------------------------*/

/**
 * \brief Return a set of edges reachable for reconstruction by the edge
 * \param[in] p_edge Edge used to determine reachable edges
 */
edgeSet DynamicNetwork::get_reachableEdgeSet(edge p_edge) const
{
	edgeSet reachableEdgeSet;
	for (auto iter = m_neighborMap.at(p_edge.first).begin(); 
		iter != m_neighborMap.at(p_edge.first).end(); ++iter)
	{
		edge e(p_edge.first,*iter);
		order(e);
		if (m_remainingEdgeMap.at(e) > 0)
		{
			reachableEdgeSet.insert(e);
		}
	}
	for (auto iter = m_neighborMap.at(p_edge.second).begin(); 
		iter != m_neighborMap.at(p_edge.second).end(); ++iter)
	{
		edge e(p_edge.second,*iter);
		order(e);
		if (m_remainingEdgeMap.at(e) > 0)
		{
			reachableEdgeSet.insert(e);
		}
	}

	return reachableEdgeSet;
}

}//end of namespace DynNet