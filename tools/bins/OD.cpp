/********************************************************/
/*                Laurent Hébert-Dufresne               */
/* Age of edges or nodes by OD of an undirected network */
/* 			     based on arXiv:1510.08542              */
/*                   Santa Fe Institute                 */
/********************************************************/
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include <time.h>
#include <set>
//#include <boost/filesystem/operations.hpp>
//#include <boost/filesystem/fstream.hpp>
//#include <boost/filesystem/path.hpp>
using namespace std;

//COMPARISON (to order set of [node tag,degree] object by degree
struct classcomp {
  bool operator() (const pair<int,int>& lhs, const pair<int,int>& rhs) const
  {return lhs.second<rhs.second;}
};

//COMPARISON (to order set of [node tag,degree] object by degree
struct classcomp2 {
  bool operator() (const pair< pair<int,int> ,int>& lhs, const pair< pair<int,int> ,int>& rhs) const
  {return lhs.second<rhs.second;}
};

//GLOBAL VARIABLES DEFINITION
vector< vector<int> > adjmat; //[neighbours]
vector<int> degree; //[degree]
multiset< pair<int,int>, classcomp> tags; //[tag,degree]
multiset< pair<int,int> > edges; //[node1, node2]
multiset< pair< pair<int,int> ,int>, classcomp2> ranked_edges; //[tag,degree]

//MAIN
int main(int argc, const char *argv[])
{
  //NETWORK
  std::string name = argv[1]; //path to edgelist
  int nodes = 0;
  if(argc>2) nodes = atoi(argv[2]); //rank output for: 0==edges, 1==nodes

//INPUT AND CONSTRUCTION OF STRUCTURES
//PREPARES INPUT & OUTPUT
std::ifstream input(name.c_str());
string line;
ifstream input0(name.c_str());
int temp = 0;
int MAX = 0;
int node1;
int node2;
if (!input0.is_open()) {
	cerr<<"error in filenames of input directory!"<<endl;;
 	return EXIT_FAILURE;
} //end if
else {
	while ( input0 >> temp ) if(temp > MAX) MAX = temp;
}
MAX = MAX+1;
int link=0;
//cout << "Number of nodes: " << MAX << endl;
adjmat.resize(MAX);
degree.resize(MAX);
if (!input.is_open()) {
	cerr<<"error in filenames of input directory!"<<endl;;
 	return EXIT_FAILURE;
} //end if
else {
    string line_buffer;
    while (getline(input, line_buffer))
    {
      	stringstream ls(line_buffer);
      	ls >> node1;
      	ls >> node2;
	// while ( input >> node1 >> std::ws >> node2 >> std::ws >> tmp >> std::ws ) {
		//if(find(adjmat[node1].begin(),adjmat[node1].end(),node2)==adjmat[node1].end() && node1 != node2) {
		adjmat[node1].push_back(node2);
		adjmat[node2].push_back(node1);
		edges.insert(make_pair(node1,node2)); 
		++degree[node1];
		++degree[node2];
		link = link+2;
		//}
	} //end while
} //end else
input.close();
for(int el=0; el<MAX; ++el) tags.insert(make_pair(el,degree[el]));

vector<int> realindegree = degree;

//ALGORITHM
vector<int> core(MAX,-1);
vector<int> pass(MAX,-1);
vector<int> age(MAX,-1);
multiset< pair<int,int>, classcomp>::iterator v;
multiset< pair<int,int>, classcomp>::iterator it1;
multiset< pair<int,int>, classcomp>::iterator it2;
pair<multiset< pair<int,int>, classcomp>::iterator,multiset< pair<int,int>, classcomp>::iterator> ret_intern;
pair<multiset< pair<int,int>, classcomp>::iterator,multiset< pair<int,int>, classcomp>::iterator> ret_pass;
int source;
int newguy;
int neighbour;
int currentpass=1;
v = tags.begin();
int maxcore = 0; int t=0;
while(tags.size()>0) {
	newguy = (*v).first;
	if(degree[newguy]>maxcore) maxcore = degree[newguy];//keeping track of the number of cores
	ret_pass = tags.equal_range(make_pair(0,degree[newguy]));//all nodes in the next layer
	int number = tags.count(make_pair(0,degree[newguy]));//number of nodes in next layer
	it1=ret_pass.first;
	for(int tmp=0; tmp<number; ++tmp,++it1) {//loops around everyone in the same layer
		source = (*it1).first; //node in the layer
		core[source] = degree[source]; //core for that node
		pass[source] = currentpass; //layer for that node
		age[source] = t;
		for(int u=0; u<adjmat[source].size(); ++u) {
			neighbour = adjmat[source][u];
			if(degree[neighbour] > degree[source]) {
				ret_intern = tags.equal_range(make_pair(neighbour,degree[neighbour]));
				for(it2=ret_intern.first; it2!=ret_intern.second; ++it2) {
					if((*it2).first==neighbour) {
						tags.erase(it2);
						break;
					}
				}
				--degree[neighbour];
				tags.insert(make_pair(neighbour,degree[neighbour]));
			}
		}
	}
	t+=number;
	it1=ret_pass.first;
	for(int tmp=0; tmp<number; ++tmp,++it1) tags.erase(it1);
	v = tags.begin(); ++currentpass;
}

if(nodes==0) {//output edges

	for(multiset< pair<int,int> >::iterator e = edges.begin(); e!=edges.end(); ++e) {
		int age0=t-age[(*e).first];
		if(t-age[(*e).second]>age0) age0=t-age[(*e).second];
		ranked_edges.insert(make_pair((*e),age0));
	}

	int time = 0;
 	int	class_id = 0;
 	int class_base_time = 0;
	multiset< pair<int,int> > class_content;

	for(multiset< pair< pair<int,int>,int> >::iterator e = ranked_edges.begin(); e!=ranked_edges.end(); ++e) {
		if (e->second != class_id)
		{
			pair<int,int> last_out_edge = make_pair(class_content.begin()->first, class_content.begin()->second);
			int out_edge_tag = -1;
			for (multiset< pair<int,int> >::iterator e2 = class_content.begin(); e2 != class_content.end(); ++e2)
			{
				if (e2->first == last_out_edge.first && e2->second == last_out_edge.second)
				{
					++out_edge_tag;
				}
				else
				{
					out_edge_tag = 0;
					last_out_edge = make_pair(e2->first, e2->second);
				}
				cout << e2->first << " " << e2->second << " " << out_edge_tag << " " << class_base_time + float(class_content.size() + 1)/2 -1 << "\n";
			}
			class_content.clear();
			class_id = (*e).second;
			class_base_time = time;
		}
		class_content.insert(e->first);
		++time;
	}
	pair<int,int> last_out_edge = make_pair(class_content.begin()->first, class_content.begin()->second);
	int out_edge_tag = -1;
	for (multiset< pair<int,int> >::iterator e2 = class_content.begin(); e2 != class_content.end(); ++e2)
	{
		if (e2->first == last_out_edge.first && e2->second == last_out_edge.second)
		{
			++out_edge_tag;
		}
		else
		{
			out_edge_tag = 0;
			last_out_edge = make_pair(e2->first, e2->second);
		}
		cout << e2->first << " " << e2->second << " " << out_edge_tag << " " << class_base_time + float(class_content.size() + 1)/2 -1 << "\n";
	}
}
else {//output nodes
	
	for(int dummy=0; dummy<core.size(); ++dummy) cout << dummy << " " << pass[dummy] << "\n";

}


return 0;
} //end main
