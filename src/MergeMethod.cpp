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
// DOI: https://...
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



#include "MergeMethod.h"

MergeMethod::MergeMethod():biLouvainMethodMurataPN(){}
MergeMethod::~MergeMethod(){}

std::map<int,std::vector<int>> MergeMethod::mergeMethodCalculation(Graph &g)
{
	std::map<int,std::vector<int>> communitiesToMerge;
	int communitiesKey[_numberCommunities];
	int startFor = 0;
	int endFor = g._lastIdPartitionV1+1;
	int maxIntersection = 0;
	int key = 0;
	std::vector<int> a;
	std::vector<int> b;

	initialCommunityDefinition(g);  						//each node belongs to its individual community
	initialCommunityNeighborsDefinition(g);					//defining neighbor communities
	double totalModularity = CoClusterMateDefinitionAllCommunities(g); 	//assign the cocluster mate community(ies) to each community

	for(int i=0;i<_numberCommunities;i++)
	{
		if(i > g._lastIdPartitionV1)
		{
			startFor = g._lastIdPartitionV1+1;
			endFor = _numberCommunities;
		}

		maxIntersection = 0;
		key = i;
		//StringSplitter::printVector(_communities[i].getCorrespondenceCommunityId());
		a =  _communities[i].getCoClusterMateCommunityId();
		std::unordered_set<int> c(a.begin(),a.end());
		//for(int j=startFor;j<endFor;j++)
		for(int j=startFor;j<i;j++)
		{
			//if(i!=j)
			//{
				b = _communities[j].getCoClusterMateCommunityId();
				int intersection = std::count_if(b.begin(),b.end(),[&](int k){return c.find(k) != c.end();});
				//set_intersection(_communities[j].getCoClusterMateCommunityId().begin(),_communities[j].getCoClusterMateCommunityId().end(),_communities[i].getCoClusterMateCommunityId().begin(),_communities[i].getCoClusterMateCommunityId().end(),back_inserter(intersection));
				if(intersection > maxIntersection)
				{
					key = j;
					maxIntersection = intersection;
				}
			//}
			b.clear();
		}
		a.clear();
		//if(i != communitiesKey[key])
		//{
			communitiesKey[i] = key;
			communitiesToMerge[key].push_back(i);
		/*}
		else
		{
			communitiesToMerge[i].push_back(i);
			communitiesKey[i] = i;
		}*/
		//std::cout << "\nCommunity: " << i << "   Key:  " << key << std::endl;
	}
	return communitiesToMerge;
}


std::map<int,std::vector<int>> MergeMethod::mergeMethodCalculationWithUpdates(Graph &g)
{
	std::map<int,std::vector<int>> communitiesToMerge;
	int startFor = 0;
	int maxIntersection = 0;
	int key = 0;
	std::vector<int> a;
	std::vector<int> b;

	initialCommunityDefinition(g);  						//each node belongs to its individual community
	initialCommunityNeighborsDefinition(g);						//defining neighbor communities
	double totalModularity = CoClusterMateDefinitionAllCommunities(g); 		//assign the cocluster community(ies) to each community
	
	for(int i=0;i<_numberCommunities;i++)
	{
		if(i > g._lastIdPartitionV1)
			startFor = g._lastIdPartitionV1+1;
		maxIntersection = 0;
		key = i;
		a =  _communities[i].getCoClusterMateCommunityId();
		std::unordered_set<int> c(a.begin(),a.end());
		for(int j=startFor;j<i;j++)
		{
			//std::cout << "\nCommunity: " << j << "\t Number Nodes: " << _communities[j].getNumberNodes() << std::endl;
			if(_communities[j].getNumberNodes()>0)
			{
				b = _communities[j].getCoClusterMateCommunityId();
				int intersection = std::count_if(b.begin(),b.end(),[&](int k){return c.find(k) != c.end();});
				//set_intersection(_communities[j].getCoClusterMateCommunityId().begin(),_communities[j].getCoClusterMateCommunityId().end(),_communities[i].getCoClusterMateCommunityId().begin(),_communities[i].getCoClusterMateCommunityId().end(),back_inserter(intersection));
				if(intersection > maxIntersection)
				{
					key = j;
					maxIntersection = intersection;
				}
				b.clear();
				//std::cout << "\nCommunity Comparison: " << j << "\t Inter: " << intersection << std::endl;
			}
		}
		a.clear();
		if(key != i) //we have communities to merge
		{
			//update communities
			updateNodeCommunity(g,g._graph[i].getId(),i,key);
			//update merged community neighbors
			updateNeighborCommunities(g,g._graph[i].getId(),i,key);
			//update merged community idCoClusters
			CoClusterMateDefinitionIDCommunity(g,key);
		}
		//std::cout << "\nCommunity: " << i << "   Key:  " << key << std::endl;
	}
	for(int i=0;i<_numberCommunities;i++)
	{
		if(_communities[i].getNumberNodes()>0)
		{
			for(unsigned int j=0;j<_communities[i].getNodes().size();j++)
				communitiesToMerge[i].push_back(_communities[i].getNodes()[j]);
		}
	}
	fromCommunitiesToNodes(g);
	//for(int i=0;i<g._numberNodes;i++)
        //        std::cout << g._graph[i].getId()<< "\t" << g._graph[i].getDegreeNode()<<std::endl;
	return communitiesToMerge;
}



void MergeMethod::initialCommunityDefinitionProvidedFile(Graph &g,const std::string &initialCommunitiesFileName)
{
        std::ifstream initialCommunitiesFile(initialCommunitiesFileName.c_str());
        std::tr1::unordered_map<int,double> nodesInCommunity;
        int numberCommunities = 0;
        int numberNodes = 0;
        int id = 0;
        std::string line;
        std::string* nodes;
        if(initialCommunitiesFile.is_open())
        {
		while(initialCommunitiesFile.good())
                {
                        getline(initialCommunitiesFile,line);
                        if(initialCommunitiesFile.eof())break;
                        if(line.find(",")!= std::string::npos)
                                nodes = StringSplitter::split(line,",",numberNodes);
                        else
                                nodes = StringSplitter::split(line,"\n",numberNodes);                        
			for(int i=0;i<numberNodes;i++)
                        {
                                id = stoi(nodes[i]);
                                g._graph[id].setCommunityId(numberCommunities);
				nodesInCommunity[g._graph[id].getId()] = g._graph[id].getDegreeNode();
                        }
                        Community community(numberCommunities,g._graph[id].getType(),nodesInCommunity);
                        _communities.push_back(community);
                        nodesInCommunity.clear();
                        numberCommunities++;
                }
                _numberCommunities = numberCommunities;
                delete[] nodes;
                initialCommunitiesFile.close();
		fromCommunitiesToNodes(g);
	        //for(int i=0;i<g._numberNodes;i++)
                //      std::cout << g._graph[i].getId()<< "\t" << g._graph[i].getDegreeNode()<<std::endl;
        }
        else
        {
                printf("\nInitial Communities File not found\n");
                exit(EXIT_FAILURE);
        }
}


void MergeMethod::mergeMethodFile(Graph &g,const std::string &inputFileName)
{
	std::stringstream initialCommunities;
	std::map<int,std::vector<int>> communitiesToMerge;
	int pos = inputFileName.find(".");
	std::string outputFileName = inputFileName.substr(0,pos) + "_InitialCommunities.txt";
	std::ifstream inputFile(outputFileName.c_str());
	if(inputFile.is_open())
		initialCommunityDefinitionProvidedFile(g,outputFileName);
	else
	{
		communitiesToMerge = mergeMethodCalculationWithUpdates(g);
        	std::ofstream outputFile;
	        outputFile.open(outputFileName.c_str(),std::ios::out|std::ios::trunc);
		for(auto it=communitiesToMerge.begin();it!=communitiesToMerge.end();++it)
		{
			for(unsigned int j=0;j<it->second.size();j++)
				initialCommunities << it->second[j] << ",";
			initialCommunities.seekp(initialCommunities.str().length()-1);
			initialCommunities << "\n";
			outputFile << initialCommunities.str();
			initialCommunities.str("");
		}
		outputFile.close();
	}
}
