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


#include "Timer.h"

double mytimer(void)
{
  static long int start = 0L, startu;

  const double million = 1000000.0;

  timeval tp;

  if (start == 0L) {
    gettimeofday(&tp, NULL);

    start = tp.tv_sec;
    startu = tp.tv_usec;
  }

  gettimeofday(&tp, NULL);

  return (static_cast<double>(tp.tv_sec - start) + (static_cast<double>(tp.tv_usec - startu) /
						    million));
}


std::string timeConverter(double time)
{
	std::stringstream timeConverted;
	long microseconds   = (long) time % 1000;
	long milliseconds   = (long) (time / 1000) % 1000;
	long seconds    = (((long) (time / 1000) - milliseconds)/1000)%60 ;
	long minutes    = (((((long) (time / 1000) - milliseconds)/1000) - seconds)/60) %60;
	long hours      = ((((((long) (time / 1000) - milliseconds)/1000) - seconds)/60) - minutes)/60;
	timeConverted << hours << " H :: " << minutes << " Min :: " << seconds << " Sec :: "<< milliseconds << " Milli :: " << microseconds << " Micro";
	return timeConverted.str();
}

