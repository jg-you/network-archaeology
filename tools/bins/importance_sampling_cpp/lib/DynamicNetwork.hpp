/**
* \file DynamicNetwork.hpp
* \brief Header file for class DynamicNetwork
* \author Guillaume St-Onge
* \version 1.0
* \date 08/11/2017
*/

#include <utility>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <boost/functional/hash.hpp>

#ifndef DYNAMICNETWORK_HPP_
#define DYNAMICNETWORK_HPP_

namespace DynNet
{//start of namespace DynNet

typedef unsigned int node;
typedef std::pair<node, node> edge;
typedef std::unordered_map<edge, unsigned int, 
	boost::hash<edge> > edgeIntMap;
typedef std::unordered_set<edge, boost::hash<edge> > edgeSet;

/**
* \class DynamicNetwork DynamicNetwork.hpp
* \brief Structure of a network in reconstruction.
*/

class DynamicNetwork
{
public:
    //Constructors
    DynamicNetwork(std::unordered_map<node,std::set<node> >& p_neighborMap,
    	edgeIntMap& p_multiplicityMap);

    //Accessors
    const std::vector<edge>& get_currentEdgeList() const
    	{return m_currentEdgeList;}
    const edgeSet& get_currentEdgeSet() const
    	{return m_currentEdgeSet;}
    const edgeSet& get_reachableEdgeSet() const
    	{return m_reachableEdgeSet;}
    edgeSet get_reachableEdgeSet(edge p_edge) const; //reachable by p_edge
    const std::unordered_map<node, unsigned int>& get_currentDegreeMap() const
    	{return m_currentDegreeMap;}
    const edgeIntMap& get_remainingEdgeMap() const
    	{return m_remainingEdgeMap;}
    double get_logweight() const
    	{return m_logweight;}
    unsigned int get_reachableEdge() const
        {return m_reachableEdge;}

    //Mutators
    void update_reachableEdge(edge p_edge);
    void init(edge p_seedEdge, double p_seedInvWeight);
    void init(edge p_seedEdge);
    void add(edge p_edge, double p_weightNorm);

private:
    //Members
    std::unordered_map<node,std::set<node> > m_neighborMap;
    std::vector<edge> m_currentEdgeList;
	edgeSet m_currentEdgeSet;
	edgeSet m_reachableEdgeSet;
	std::unordered_map<node, unsigned int> m_currentDegreeMap;
	edgeIntMap m_remainingEdgeMap;
	double m_logweight;
    unsigned int m_reachableEdge;
};

}//end of namespace DynNet

#endif /* DYNAMICNETWORK_HPP_ */