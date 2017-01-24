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


#include "biLouvainMethodMurataPN.h"

biLouvainMethodMurataPN::biLouvainMethodMurataPN():biLouvainMethod(){}
biLouvainMethodMurataPN::~biLouvainMethodMurataPN(){}

double biLouvainMethodMurataPN::murataModularityArgMax(Graph &g,int &communityId, int possibleCoClusterMateId)
{
	double al = (double)_communities[communityId].getDegreeCommunity();
	double am = (double)_communities[possibleCoClusterMateId].getDegreeCommunity();
	double murataModularity = 0.0;
	//printf("\nCommunity ID: %d \n",communityId);
	//printf("AL: %f AM: %f  Edges V1: %d Edges V2: %d  Total edges: %d\n",al,am,_number_edges_V1,_number_edges_V2,_number_edges);
	if(_communities[communityId].getDescription()=="V1")
	{	al = al/(2*g._weightEdgesV1);
		am = am /(2*g._weightEdgesV2);
	}
	else
	{	al = al /(2*g._weightEdgesV2);
		am = am / (2*g._weightEdgesV1);
	}

	double elm = ((double)1/(2*g._weightEdges)) * calculateEdgesBetweenCommunitiesMap(g,communityId,possibleCoClusterMateId);
	//printf("AM: %f  AL: %f  Edges V1: %d Edges V2: %d  Total edges: %d Edges between Communities:%f \n",am,al,_number_edges_V1,_number_edges_V2,_number_edges,calculateEdgesBetweenCommunitiesMap(communityId,possibleCoClusterMateId));
	murataModularity = elm-(al*am);
	//printf("Possible cocluster mate ID: %d  Murata Modularity: %f \n",possibleCoClusterMateId,murataModularity);
	return murataModularity;
}

double biLouvainMethodMurataPN::murataModularityWithChanges(Graph &g,MetaNode &node, int &communityId, int possibleCoClusterMateId,int &newCommunityId,int &option)
{
	double al = 0.0;
	double am = 0.0;
	double elm = 0.0;
	double murataModularity = 0.0;
	int nodeCommunity = -1;

	switch(option)
	{
		case 1: // For Ci community
		{
			al = (double)_communities[communityId].getDegreeCommunityWithoutNode(node.getId());
			am = (double)_communities[possibleCoClusterMateId].getDegreeCommunity();
			elm = ((double)1/(2*g._weightEdges)) * (calculateEdgesBetweenCommunitiesMap(g,communityId,possibleCoClusterMateId)-node.getWeightEdgesToNeighborCommunity(possibleCoClusterMateId));
			break;
		}
		case 2: // For Cj community
		{
			al = (double)_communities[communityId].getDegreeCommunity() + node.getDegreeNode();
			am = (double)_communities[possibleCoClusterMateId].getDegreeCommunity();
			elm = ((double)1/(2*g._weightEdges)) * (calculateEdgesBetweenCommunitiesMap(g,communityId,possibleCoClusterMateId)+node.getWeightEdgesToNeighborCommunity(possibleCoClusterMateId));
			break;
		}
		case 3: case 4: // For Di community which is a cocluster of Ci
		{
			nodeCommunity = node.getCommunityId();
			al = (double)_communities[communityId].getDegreeCommunity();
			if(possibleCoClusterMateId == nodeCommunity)
			{
				am = (double)_communities[possibleCoClusterMateId].getDegreeCommunity() - (double)node.getDegreeNode();
				elm = ((double)1/(2*g._weightEdges)) * (calculateEdgesBetweenCommunitiesMap(g,communityId,possibleCoClusterMateId)-node.getWeightEdgesToNeighborCommunity(communityId));
			}
			else if(possibleCoClusterMateId == newCommunityId)
			{
				am = (double)_communities[possibleCoClusterMateId].getDegreeCommunity() + (double)node.getDegreeNode();
				elm = ((double)1/(2*g._weightEdges)) * (calculateEdgesBetweenCommunitiesMap(g,communityId,possibleCoClusterMateId)+node.getWeightEdgesToNeighborCommunity(communityId));
			}
			else
			{
				am = (double)_communities[possibleCoClusterMateId].getDegreeCommunity();
				elm = ((double)1/(2*g._weightEdges)) * (calculateEdgesBetweenCommunitiesMap(g,communityId,possibleCoClusterMateId));
			}
			break;
		}
	}
	if(_communities[communityId].getDescription()=="V1")
	{	al = al/(2*g._weightEdgesV1);
		am = am /(2*g._weightEdgesV2);
	}
	else
	{	al = al /(2*g._weightEdgesV2);
		am = am / (2*g._weightEdgesV1);
	}
	//printf("AM: %f  AL: %f  Total edges: %d Edges between Communities:%d \n",am,al,_number_edges,calculateEdgesBetweenCommunities(g,communityId,possibleCoClusterMateId));
	murataModularity = elm-(al*am);
	//printf("Possible cocluster mate ID: %d  Murata Modularity: %f \n",possibleCoClusterMateId,murataModularity);
	return murataModularity;
}

newDataCommunityVector biLouvainMethodMurataPN::murataCalculationCoClusterMates(Graph &g,int communityId,const std::vector<int> &possibleCoClusterMates)
{
	double murataModularity = 0.0;
        double maxMurataModularity = 0.0;
	std::vector<int> coClusterMateCommunityId;
	newDataCommunityVector result;

	for(unsigned int j=0;j<possibleCoClusterMates.size();j++)
        {
	        murataModularity = murataModularityArgMax(g,communityId,possibleCoClusterMates[j]);
                //printf("Possible cocluster mate ID: %d  Murata Modularity: %f \n",possibleCoClusterMates[j],murataModularity);
                if(j==0)
              	{
                	maxMurataModularity = murataModularity;
                	coClusterMateCommunityId.push_back(possibleCoClusterMates[0]);
                }
                else
                {
                      	if(murataModularity > maxMurataModularity)
                     	{
                            	maxMurataModularity = murataModularity;
                              	coClusterMateCommunityId.clear();
                            	coClusterMateCommunityId.push_back(possibleCoClusterMates[j]);
                        }
                        else if(murataModularity == maxMurataModularity)
                                coClusterMateCommunityId.push_back(possibleCoClusterMates[j]);
                }
                //printf("From: %d  To: %d  AL: %f  AM: %f  ELM: %f  Murata M: %f   Max Murata M: %f \n",_communities[i]etId(),possibleCoClusterMates[i],al,am,elm,murataModularity,maxMurataModularity);
        }
	result.coClusterMateCommunityId = coClusterMateCommunityId;
	result.newModularityContribution = maxMurataModularity;
	coClusterMateCommunityId.clear();
	return result;
}


double biLouvainMethodMurataPN::CoClusterMateDefinitionAllCommunities(Graph &g, int start, int end)
{
	std::vector<int> possibleCoClusterMates;
	std::vector<int> coClusterMateCommunityId;
	double totalPartitioningModularityCalculated = 0.0;
	double maxMurataModularity = 0.0;	
	newDataCommunityVector communityModularity;

	for(int i=start;i<end;i++)
	{
		if(_communities[i].getNumberNodes()>0)
		{
			//Part I: Find the possible cocluster mates
			possibleCoClusterMates = findNeighborCommunitiesMap(g,_communities[i].getId());

			//Part II: Apply Murata+ calculation
			communityModularity = murataCalculationCoClusterMates(g,_communities[i].getId(),possibleCoClusterMates);
			coClusterMateCommunityId = communityModularity.coClusterMateCommunityId;
			maxMurataModularity = communityModularity.newModularityContribution;
			//printf("\n Comm: %d  Max Mod: %f \n",i,maxMurataModularity);
			_communities[i].setModularityContribution(maxMurataModularity);
			totalPartitioningModularityCalculated += maxMurataModularity;

			//Part III: Assign the collection of possible cocluster mates to the community and add the contribution of each community to the total modularity
			//-2 Empty community
			//-1 Nodes inside the community don't have neighbors
			if(coClusterMateCommunityId.size()>0)
				_communities[i].setCoClusterMateCommunityId(coClusterMateCommunityId);
			else
			{
				coClusterMateCommunityId.push_back(-1);
				_communities[i].setCoClusterMateCommunityId(coClusterMateCommunityId);
			}
			possibleCoClusterMates.clear();
			coClusterMateCommunityId.clear();
		}
		else
		{
			coClusterMateCommunityId.push_back(-2);
			_communities[i].setCoClusterMateCommunityId(coClusterMateCommunityId);
			_communities[i].setModularityContribution(0.0);
			coClusterMateCommunityId.clear();
		}
	}
	return totalPartitioningModularityCalculated;
}


void biLouvainMethodMurataPN::CoClusterMateDefinitionIDCommunity(Graph &g,int &communityId)
{
	std::vector<int> possibleCoClusterMates;
	std::vector<int> coClusterMateCommunityId;
	double maxMurataModularity = 0.0;
	newDataCommunityVector communityModularity;

	if(_communities[communityId].getNumberNodes()>0)
	{
		//Part I: Find the possible cocluster mates
		possibleCoClusterMates = findNeighborCommunitiesMap(g,communityId);

		//Part II: Apply Murata+ calculation
		communityModularity = murataCalculationCoClusterMates(g,communityId,possibleCoClusterMates);
                coClusterMateCommunityId = communityModularity.coClusterMateCommunityId;
                maxMurataModularity = communityModularity.newModularityContribution;
                //printf("\n Comm: %d  Max Mod: %f \n",i,maxMurataModularity);
		_communities[communityId].setModularityContribution(maxMurataModularity);

		//Part III: Return the collection of possible cocluster mates to the community
		if(coClusterMateCommunityId.size()>0)
			_communities[communityId].setCoClusterMateCommunityId(coClusterMateCommunityId);
		else
		{
			coClusterMateCommunityId.push_back(-1);
			_communities[communityId].setCoClusterMateCommunityId(coClusterMateCommunityId);
		}
		possibleCoClusterMates.clear();
		coClusterMateCommunityId.clear();
	}
	else
	{
		coClusterMateCommunityId.push_back(-2);
		_communities[communityId].setCoClusterMateCommunityId(coClusterMateCommunityId);
		_communities[communityId].setModularityContribution(0.0);
		coClusterMateCommunityId.clear();
	}
}


newDataCommunity biLouvainMethodMurataPN::CoClusterMateDefinitionPrecalculation(Graph &g,MetaNode &node, int &communityId, int &newCommunityId,int option)
{
	std::vector<int> possibleCoClusterMates;
	std::vector<int> coClusterMateCommunityId;
	int communityContainingNodeId = 0;
	double maxMurataModularity = 0.0;
	double murataModularity = 0.0;
	std::stringstream changes;
	newDataCommunity result;
	struct timeval t1,t2,t3,t4,t5,t6,t7,t8;

	if(_communities[communityId].getNumberNodes()>0)
	{
		//Part I: Find possible cocluster mates

		if(option==1)// For Ci community
		{
			gettimeofday(&t1,NULL);
			possibleCoClusterMates = findNeighborCommunitiesWithoutNodeMap(g,communityId,node.getId());
			gettimeofday(&t2,NULL);
			precalculationCiTime += (t2.tv_sec - t1.tv_sec)*1000000 + (t2.tv_usec - t1.tv_usec);
		}
		else if(option==2)// For Cj community
		{
			gettimeofday(&t3,NULL);
			possibleCoClusterMates = findNeighborCommunitiesMap(g,communityId);
			std::vector<int> temp = node.getNeighborCommunities();
			possibleCoClusterMates.insert(possibleCoClusterMates.end(),temp.begin(),temp.end());
			sort(possibleCoClusterMates.begin(),possibleCoClusterMates.end());
			possibleCoClusterMates.erase(unique(possibleCoClusterMates.begin(),possibleCoClusterMates.end()),possibleCoClusterMates.end());
			gettimeofday(&t4,NULL);
			precalculationCjTime += (t4.tv_sec - t3.tv_sec)*1000000 + (t4.tv_usec - t3.tv_usec);
		}
		else if((option==3)||(option==4))// For Di community which is a cocluster of Ci
		{
			gettimeofday(&t5,NULL);
			std::vector<int> temp;
			for(int i=0;i<_communities[communityId].getNumberNodes();i++)
			{
				for(int j=0;j<g._graph[_communities[communityId].getNodes()[i]].getNumberNeighbors();j++)
				{
					if(g._graph[_communities[communityId].getNodes()[i]].getNeighbors()[j]==node.getId()) communityContainingNodeId = newCommunityId;
					else communityContainingNodeId = g._graph[g._graph[_communities[communityId].getNodes()[i]].getNeighbors()[j]].getCommunityId();
					std::vector<int>::iterator position = find(temp.begin(),temp.end(),communityContainingNodeId);
					if (position == temp.end())
						temp.push_back(communityContainingNodeId);
				}
				possibleCoClusterMates.insert(possibleCoClusterMates.end(),temp.begin(),temp.end());
				temp.clear();
			}
			sort(possibleCoClusterMates.begin(),possibleCoClusterMates.end());
			possibleCoClusterMates.erase(unique(possibleCoClusterMates.begin(),possibleCoClusterMates.end()),possibleCoClusterMates.end());
			gettimeofday(&t6,NULL);
			precalculationDTime += (t6.tv_sec - t5.tv_sec)*1000000 + (t6.tv_usec - t5.tv_usec);
		}

		//Part II: Apply Murata+ calculation
		gettimeofday(&t7,NULL);
		for(unsigned int j=0;j<possibleCoClusterMates.size();j++)
        	{
                	murataModularity = murataModularityWithChanges(g,node,communityId,possibleCoClusterMates[j],newCommunityId,option);
	                //printf("Possible cocluster mate ID: %d  Murata Modularity: %f \n",possibleCoClusterMates[j],murataModularity);
        	        if(j==0)
                	{
                        	maxMurataModularity = murataModularity;
	                        coClusterMateCommunityId.push_back(possibleCoClusterMates[0]);
        	        }
                	else
	                {	
        	                if(murataModularity > maxMurataModularity)
                	        {
                        	        maxMurataModularity = murataModularity;
                                	coClusterMateCommunityId.clear();
	                                coClusterMateCommunityId.push_back(possibleCoClusterMates[j]);
        	                }
                	        else if(murataModularity == maxMurataModularity)
                        	        coClusterMateCommunityId.push_back(possibleCoClusterMates[j]);
               		 }
	                //printf("From: %d  To: %d  AL: %f  AM: %f  ELM: %f  Murata M: %f   Max Murata M: %f \n",_communities[i]etId(),possibleCoClusterMates[i],al,am,elm,murataModularity,maxMurataModularity);
       		 }
		//printf("\n Max Mod: %f \n",maxMurataModularity);

		//Part III: Return the collection of possible cocluster mates to the community
		//-2 Empty community
		//-1 Nodes inside the community don't have neighbors
		if(coClusterMateCommunityId.size()==0)
		{
			if(_communities[communityId].getNodesWithoutNode(node.getId()).size()>0)
				coClusterMateCommunityId.push_back(-1);
			else
				coClusterMateCommunityId.push_back(-2);
		}
		gettimeofday(&t8,NULL);
		premurataTime += (t8.tv_sec - t7.tv_sec)*1000000 + (t8.tv_usec - t7.tv_usec);
	}
	else
		coClusterMateCommunityId.push_back(-2);
	changes << communityId << ":";
	for(unsigned int i=0;i<coClusterMateCommunityId.size();i++)
		changes << coClusterMateCommunityId[i] << ",";
	changes.str(changes.str().substr(0,changes.str().length()-1));
	result.coClusterMateCommunityId = changes.str();
	result.newModularityContribution = maxMurataModularity;
	//printf("Max:%f  Back:%f\n",maxMurataModularity, =.at(coClusterMateCommunityId.size()-1));
	return result;
}
