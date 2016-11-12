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

/*
# Community.h
# Represents a set of nodes of the same type that are strongly connected through the strength of
# their shared connections to nodes of the other type.
*/


#ifndef COMMUNITY_H_
#define COMMUNITY_H_


#include "Header.h"

struct newDataCommunity
{
	std::string coClusterMateCommunityId;
	double newModularityContribution;
};

struct newDataCommunityVector
{
        std::vector<int> coClusterMateCommunityId;
        double newModularityContribution;
};

class Community
{
	private:
		int _id;
		std::string _description;
		double _modularityContribution;
		std::vector<int> _coClusterMateCommunityId;
		std::tr1::unordered_map<int,double> _nodes;

	public:
		Community();
		Community(int id, std::string description, std::tr1::unordered_map<int,double> nodes);

		/* Get functions */
		int getId();
		std::string getDescription();
		double getModularityContribution();
		std::vector<int> getNodes();
		std::vector<int> getCoClusterMateCommunityId();
		int getNumberNodes();
		double getDegreeCommunity();
		double getDegreeCommunityWithoutNode(int nodeId);
		std::vector<int> getNodesWithoutNode(int nodeId);

		/* Set procedures */
		void setId(int id);
		void setDescription(std::string description);
		void setCoClusterMateCommunityId(std::vector<int> coClusterCommunityId);
		void setNodes(std::tr1::unordered_map<int,double> nodes);
		void setModularityContribution(double modularityContribution);
		void addNode(int nodeId, double nodeDegree);
		void deleteNode(int nodeId);
};


#endif /* COMMUNITY_H_ */