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
// MeshingSpheresMethods.h                                        Last modified (07/12/2018) //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _SPHERES_METHODS_H_
#define _SPHERES_METHODS_H_

#include <cmath>
#include <vector>
#include <list>
#include <utility>
#include <stack>
#include <map>
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

#include <time.h>


#include "MeshingSmartTree.h"
#include "MeshingGeometricalMethods.h"
#include "MeshingMemoryHandler.h"

class MeshingSpheresMethods
{
    
public:
    
	MeshingSpheresMethods();
    
    ~MeshingSpheresMethods();
    
	bool point_covered(double* point, MeshingSmartTree* spheres, double Lip);

	bool point_covered(double* point, double* sphere, double alpha_coverage);

	bool point_covered(double* point, MeshingSmartTree* spheres, double Lip, double alpha_coverage);

	bool edge_covered(double** edge_corners, MeshingSmartTree* spheres, double Lip, double alpha_coverage);
	
	bool face_covered(double** face_corners, MeshingSmartTree* spheres, double Lip, double alpha_coverage);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int impose_lipschitz_continuity_brute(MeshingSmartTree* spheres, double Lip, bool &shrunk);

	int impose_lipschitz_continuity(MeshingSmartTree* spheres, double Lip, bool &shrunk);

	int get_overlapping_spheres(double* sphere, MeshingSmartTree* spheres, double Lip,
		                        size_t &num_overlapping_spheres, size_t* &overlapping_spheres);

	int resolve_sliver(size_t num_spheres, double** spheres, bool* fixed, size_t &sliver_sphere_index, double &sliver_sphere_radius);

	int get_power_vertex(size_t num_dim, double* sphere_i, double* sphere_j, double* pv);

	bool get_power_vertex(size_t num_dim, size_t num_spheres, double** spheres, double* pv);


private:
	int get_normal_component(size_t num_dim, size_t num_basis, double** basis, double* vect, double &norm);

private:
	MeshingGeometricalMethods _geom;
	MeshingRandomSampler      _rsampler;
	MeshingMemoryHandler      _memo;
};

#endif

