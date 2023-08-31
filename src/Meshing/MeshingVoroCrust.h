///////////////////////////////////////////////////////////////////////////////////////////////
//                                VOROCRUST-MESHING 1.0                                      //
// Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).        //
// Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain  //
// rights in this software.                                                                  //
//                                                                                           //
// Redistribution and use in source and binary forms, with or without modification, are      //
// permitted provided that the following conditions are met:                                 //  
//                                                                                           //
// 1. Redistributions of source code must retain the above copyright notice, this list of    //
// conditions and the following disclaimer.                                                  //
//                                                                                           //
// 2. Redistributions in binary form must reproduce the above copyright notice, this list    //
// of conditions and the following disclaimer in the // documentation and/or other materials //
// provided with the distribution.                                                           //
//                                                                                           //
// 3. Neither the name of the copyright holder nor the names of its contributors may be      //
// used to endorse or promote products derived from this software without specific prior     //
// written permission.                                                                       //
//-------------------------------------------------------------------------------------------//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY       //
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF   //
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE//
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, //
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        //
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)    //
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR  //
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        //
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                              //
///////////////////////////////////////////////////////////////////////////////////////////////
//                                     Author                                                //
//                                Mohamed S. Ebeida                                          //
//                                msebeid@sandia.gov                                         //
///////////////////////////////////////////////////////////////////////////////////////////////
//  MeshingVoroCrust.h                                            Last modified (08/03/2022) //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _VOROCRUST_H_
#define _VOROCRUST_H_

#include <string>

#include <cmath>
#include <list>
#include <vector>
#include <utility>
#include <stack>
#include <map>
#include <memory>
#include <set>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <cstddef>
#include <cstdlib>
#include <limits>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <iomanip>

#include <time.h>

#include "Config.h"
#include "MeshingConsole.h"
#include "MeshingOptionParser.h"
#include "MeshingTimer.h"

#if defined USE_OPEN_MP
#include <omp.h>
#endif


#ifndef DBL_MAX
#define DBL_MAX  1.7976931348623158e+308 // max double value
#endif

#ifndef PI
#define PI  3.141592653589793
#endif

// extern std::shared_ptr<VoroCrust_OptionParser> options;

extern MeshingConsole vcm_cout;


#include "MeshingSmartTree.h"
#include "MeshingRandomSampler.h"

#include "MeshingVoroCrustSampler.h"
#include "MeshingVoronoiMesher.h"

class MeshingVoroCrust
{

public:

	MeshingVoroCrust(std::string input_filename);

	~MeshingVoroCrust();

	int execute();

private:

	int execute_vc(int num_threads, size_t num_points, double** points, size_t num_faces, size_t** faces, MeshingSmartTree* input_sizing_function,
		                double rmin, double rmax, double Lip, double smooth_angle_threshold, bool generate_vcg_file,bool generate_exodus_file,
		                bool generate_monitoring_points, bool impose_monitoring_points);

	int generate_interior_seeds(MeshingSmartTree* seeds_tree, MeshingSmartTree* surface_spheres_tree, MeshingSmartTree* edge_spheres_tree, MeshingSmartTree* corner_spheres_tree,
		                        MeshingSmartTree* sz_function_tree, int num_threads, double Lip, double rmax,
		                        bool impose_monitoring_points, bool generate_monitoring_points,
		                        size_t& num_seeds, double*& seeds, size_t*& seeds_region_id, double*& seeds_sizing);

	int generate_explicit_mesh(int num_threads, size_t num_seeds, double* seeds, size_t* seeds_region_id, double* seeds_sizing,
		                       size_t num_surface_seeds, bool generate_vcg_file, bool generate_exodus_file);

	int generate_mesh_from_seeds_file(std::string file_name, int num_threads, bool generate_vcg_file, bool generate_exodus_file);

	int read_seeds_csv(std::string file_name, size_t& num_seeds, double*& seeds, size_t*& seeds_region_id, double*& seeds_sizing);

	int get_tokens(std::string line, char separator, std::vector<std::string>& tokens);


	MeshingRandomSampler _rsampler;
	MeshingMemoryHandler _memo;

	std::string input_filename;

};

#endif
