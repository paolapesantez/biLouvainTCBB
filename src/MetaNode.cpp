// **************************************************************************************************
// biLouvain: A C++ library for bipartite graph community detection
// Paola Gabriela Pesantez-Cabrera, Ananth Kalyanaraman
//      (p.pesantezcabrera@wsu.edu, ananth@eecs.wsu.edu)
// Washington State University
//
// For citation, please cite the following paper:
// Pesantez, Paola and Kalyanaraman, Ananth, "Detecting Communities in Biological
// Bipartite Networks," Proc. ACM Conference on Bioinformatics, Computational Biology,
// and Health Informatics (ACM-BCB), Seattle, WA, October 2-5, 2016, In press,
// DOI: http://dx.doi.org/10.1145/2975167.2975177.
//
// **************************************************************************************************
// Copyright (c) 2016. Washington State University ("WSU"). All Rights Reserved.
// Permission to use, copy, modify, and distribute this software and its documentation
// for educational, research, and not-for-profit purposes, without fee, is hereby
// granted, provided that the above copyright notice, this paragraph and the following
// two paragraphs appear in all copies, modifications, and distributions. For
// commercial licensing opportunities, please contact The Office of Commercialization,
// WSU, 280/286 Lighty, PB Box 641060, Pullman, WA 99164, (509) 335-5526,
// commercialization@wsu.edu<mailto:commercialization@wsu.edu>, https://commercialization.wsu.edu/

// IN NO EVENT SHALL WSU BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
// OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
// THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF WSU HAS BEEN ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// WSU SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND
// ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". WSU HAS NO
// OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
// **************************************************************************************************


#include "MetaNode.h"

MetaNode::MetaNode(){}

MetaNode::MetaNode(int id, std::string type, std::vector<Node> nodes,std::tr1::unordered_map<int,double> neighbors,int communityId)
{
	_idGraph   = id;
	_type = type;
	_nodes = nodes;
	_neighbors = neighbors;
	_communityId = communityId;
}

/* Get functions */
int MetaNode::getId()
{
	return _idGraph;
}

std::string MetaNode::getType()
{
	return _type;
}

int MetaNode::getCommunityId()
{
	return _communityId;
}

int MetaNode::getNumberNodes()
{
	return _nodes.size();
}

std::vector<Node> MetaNode::getNodes()
{
	return _nodes;
}

std::vector<Node> MetaNode::getNodesSorted()
{
	sort(_nodes.begin(),_nodes.end(),myobject);
	return _nodes;
}

std::vector<int> MetaNode::getNeighbors()
{
	std::vector<int> result;
	for(auto it=_neighbors.begin();it!=_neighbors.end();++it)
		result.push_back(it->first);
	return result;
}

std::vector<int> MetaNode::getNeighborsSorted()
{
	std::vector<int> result;
	for(auto it=_neighbors.begin();it!=_neighbors.end();++it)
		result.push_back(it->first);
	sort(result.begin(),result.end(),myobject2);
	return result;
}

std::vector<int> MetaNode::getNeighborsWithoutNode(int nodeId)
{
	std::vector<int> result;
	for(auto it=_neighbors.begin();it!=_neighbors.end();++it)
	{
		if(it->first != nodeId)
			//cout << it->first << endl;
			result.push_back(it->first);
	}
	return result;
}

int MetaNode::getNumberNeighbors()
{
	return _neighbors.size();
}

double MetaNode::getDegreeNode()
{
	double result = 0.0;
	for(auto it=_neighbors.begin();it!=_neighbors.end();++it) 
                result += it->second;
	return result;
}

int MetaNode::getNumberNeighborCommunities()
{
	return _neighborCommunities.size();
}

double MetaNode::getDegreeNeighborCommunities()
{
	double result = 0.0;
	for(auto it=_neighborCommunities.begin();it!=_neighborCommunities.end();++it) 
		result += it->second;
	return result;
}

double MetaNode::getWeightEdgesToNeighborCommunity(int communityId)
{
	auto search = _neighborCommunities.find(communityId);
	if(search != _neighborCommunities.end())return _neighborCommunities[communityId];
	else return 0.0;
}

double MetaNode::getWeightNeighbor(int id)
{
	auto search = _neighbors.find(id);
	if(search != _neighbors.end())return _neighbors[id];
	else return 0.0;
}

std::vector<int> MetaNode::getNeighborCommunities()
{
	std::vector<int> result;
	for(auto it=_neighborCommunities.begin();it!=_neighborCommunities.end();++it)
		result.push_back(it->first);
	return result;
}


/* Set procedures */
void MetaNode::setId(int id)
{
	_idGraph = id;
}

void MetaNode::setType(std::string type)
{
	_type = type;
}

void MetaNode::setNodes(std::vector<Node> nodes)
{
	_nodes = nodes;
}

void MetaNode::setNeighbors(std::tr1::unordered_map<int,double> neighbors)
{
	_neighbors = neighbors;
}

void MetaNode::setCommunityId(int communityId)
{
	_communityId = communityId;
}

void MetaNode::setNeighborCommunities(std::tr1::unordered_map<int,double> neighborCommunities)
{
	_neighborCommunities = neighborCommunities;
}

void MetaNode::deleteNeighborCommunity(int communityId)
{
	auto search = _neighborCommunities.find(communityId);
	if(search != _neighborCommunities.end())
	{
		if(_neighborCommunities[communityId] == 1)_neighborCommunities.erase(communityId);
		else if (_neighborCommunities[communityId] > 1) _neighborCommunities[communityId]--;
	}
}

void MetaNode::deleteNeighborCommunityWeight(int communityId, double weight)
{
	auto search = _neighborCommunities.find(communityId);
	if(search != _neighborCommunities.end())
	{
		if(_neighborCommunities[communityId] == weight) _neighborCommunities.erase(communityId);
		else if (_neighborCommunities[communityId] > weight) _neighborCommunities[communityId] -=  weight;
	}
}

void MetaNode::addNeighbor(int idNeighbor,double weight)
{
	if(_neighbors.find(idNeighbor) == _neighbors.end()) _neighbors[idNeighbor] = weight;
}

void MetaNode::addNeighborCommunity(int communityId)
{
	auto search = _neighborCommunities.find(communityId);
	if(search != _neighborCommunities.end())_neighborCommunities[communityId]++;
	else _neighborCommunities[communityId] = 1;
}

void MetaNode::addNeighborCommunityWeight(int communityId, double weight)
{
	auto search = _neighborCommunities.find(communityId);
	if(search != _neighborCommunities.end())_neighborCommunities[communityId] += weight;
	else _neighborCommunities[communityId] = weight;
}

int MetaNode::findNeighborCommunity(int communityId)
{
	int result = 0;
	auto search = _neighborCommunities.find(communityId);
	if(search != _neighborCommunities.end())result = 1;
	return result;
}



