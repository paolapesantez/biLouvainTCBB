// **************************************************************************************************
// biLouvain: A C++ library for bipartite graph community detection
// Paola Gabriela Pesantez-Cabrera, Ananth Kalyanaraman  
//	(p.pesantezcabrera@wsu.edu, ananth@eecs.wsu.edu)
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

#include "Header.h"
#include "PreProcessInputBipartiteGraph.h"
#include "Graph.h"
#include "LoadGraph.h"
#include "MergeMethod.h"
#include "biLouvainMethod.h"
#include "biLouvainMethodMurataPN.h"
#include "Timer.h"


static std::string inputFileName = "";
static std::string initialCommunitiesFileName = "";
static std::string outputFileName = "";
static std::string delimiter = "\t";
static int optionOrder = 3;
static int merge = 1;
static double cutoffIterations = 0.01;
static double cutoffPhases =  0.0;
static int flag;
static void parseCommandLine(const int argc, char * const argv[]);

static struct option longopts[] = {
   { "input",		required_argument,0,'i'},
   { "delimeter",	required_argument,0,'d'},
   { "ci",		required_argument,&flag,1},
   { "cp",		required_argument,&flag,2},
   { "order",		required_argument,&flag,3},
   { "initial",		required_argument,&flag,4},
   { "merge",		required_argument,&flag,5},
   { "output",		required_argument,0,'o'},
   { 0, 0, 0, 0 }
};


int main(int argc, char *argv[])
{
	try
	{
		parseCommandLine(argc, argv);
		struct timeval startTime,endTime;	
		std::ifstream infile(inputFileName.c_str());
		//std::cout<<inputFileName<<std::endl;
		if(infile.is_open()==true)
		{
			std::string bipartiteFileName;
                        std::tr1::unordered_map<int,std::string> bipartiteOriginalEntities;
                        if(inputFileName.find("_bipartite") == std::string::npos)
                        {
				int pos = inputFileName.find(".");
                                bipartiteFileName = inputFileName.substr(0,pos)+"_bipartite.txt";
                                bipartiteOriginalEntities=PreProcessInputBipartiteGraph::preProcessingGraphData(inputFileName,delimiter);
                        }
                        else
                        {
                                bipartiteFileName = inputFileName;
                                bipartiteOriginalEntities=PreProcessInputBipartiteGraph::readDictionaryFile(inputFileName);
                        }
                        double loadGraphTime = 0.0;
                        int pass = -1;
                        gettimeofday(&startTime,NULL);
			Graph* graph;
                        pass = LoadGraph::loadBipartiteGraphFromFile(graph,bipartiteFileName);
                        gettimeofday(&endTime,NULL);
                        loadGraphTime = (endTime.tv_sec - startTime.tv_sec)*1000000 + (endTime.tv_usec - startTime.tv_usec);
                        if (pass == 0)
                         {                            
                                std::cout << "\n ::: Done Loading Bipartite Graph :::\n \n";
				MergeMethod m;
                                if(merge == 1)
                                        m.mergeMethodFile(*graph,bipartiteFileName);
                                else
                                {
                                        if(initialCommunitiesFileName.empty()==false)
                                                m.initialCommunityDefinitionProvidedFile(*graph,initialCommunitiesFileName);
                                }
				std::cout << "\n ::: Starting biLouvain Algorithm :::\n \n";
				biLouvainMethodMurataPN biLouvain;
				gettimeofday(&startTime,NULL);							
				biLouvain.biLouvainMethodAlgorithm(*graph,cutoffIterations,cutoffPhases,optionOrder,bipartiteOriginalEntities,bipartiteFileName,outputFileName);
				gettimeofday(&endTime,NULL);	
				double biLouvainAlgorithmTime = (endTime.tv_sec - startTime.tv_sec)*1000000 + (endTime.tv_usec - startTime.tv_usec);
				biLouvain.printTimes(biLouvainAlgorithmTime,loadGraphTime);						
				graph->destroyGraph();
			}
			else
			{
				printf("\n There was a problem reading the graph input file\n");
				exit(EXIT_FAILURE);
			}
		}
		else
                {
                        printf("\n Input file was not found\n");
                        exit(EXIT_FAILURE);
                }
	}
	catch (...)
	{
		std::cout << "Unknown error Main" << std::endl;
	}
	std::cout << std::endl << "::: biLouvain Method has finished :::" << std::endl;
	return 0;
}

void printUsage()
{
	 printf("Usage: -i {inputFile} -d {delimiter (\" \",\",\",\"\\t\")} [-ci {cutoff iterations(default=0.01)} -cp {cutoff phases(default=0.0)} -order {1:Sequential, 2:Alternate, 3:Random(default=3)} -initial {initialCommuitiesFile}(default=\"\") -merge {0/1 flag(default=1)} -o {outputFileName(default=input_Results*)}]\n");  
         exit(EXIT_FAILURE);
}


void parseCommandLine(const int argc, char * const argv[])
{
	int c=0,indexPtr=0,prevInd;
	while(prevInd=optind,(c = getopt_long_only(argc, argv, "i:d:o:", longopts, &indexPtr)) != -1) {
		//printf("%d \t %d \t %c\n",prevInd,optind,c);
		if((optind==prevInd+2) && (*optarg=='-'))
		{
			printf("You forgot to provide an argument %s\n");
			printUsage();
		}
		switch (c) {
		    case 'i':
		        inputFileName = optarg;
			if(inputFileName.empty())
			{	printf("No input file name provided \n");
				printUsage();
			}
		    	break;
		    case 'd':
		        delimiter = optarg;
			if(delimiter.empty())
			{	printf("No delimiter provided \n");	
				printUsage();
			}
			else if((delimiter != " ")&&(delimiter != "\\t")&&(delimiter != ","))
			{
				printf("Unknown delimiter provided \n");
		                printUsage();		
			}
			else if(delimiter == "\\t")
				delimiter = "\t";
			break;
		    case 'o':
			if(optarg != NULL)
				outputFileName = optarg;
			break;	
		    case 0:
			if(*(longopts[indexPtr].flag)==1)
			{
				if(optarg != NULL)
                        	      	cutoffIterations = atof(optarg);
			}
			else if(*(longopts[indexPtr].flag)==2)
                        {
                                if(optarg != NULL)
                                        cutoffPhases = atof(optarg);
                        }
			else if(*(longopts[indexPtr].flag)==3)
                        {
                                if(optarg != NULL)
                                        optionOrder = atoi(optarg);
                        }
			else if(*(longopts[indexPtr].flag)==4)
                        {
                                if(optarg != NULL)
                                        initialCommunitiesFileName = optarg;
                        }
			else if(*(longopts[indexPtr].flag)==5)
                        {
                                if(optarg != NULL)
                                        merge = atoi(optarg);
                        }
			break;
		    case ':':
			printUsage;
			break;
		    case '?':
			/*if((optopt == 'i') || (optopt == 'd'))
				printf("Argument is mandatory for -%c\n",optopt);
			else if(isprint(optopt))
				printf("You have provided an unknown option -%c\n",optopt);
			else
				printf("Unknown Error-0x%08x\n",optopt);			
			*/
			printUsage();
			break;
		    default:   
			exit(0);
		}
	}
}
