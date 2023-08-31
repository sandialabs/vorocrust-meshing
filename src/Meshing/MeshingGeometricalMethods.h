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
//  MeshingGeometricalMethods.h                                   Last modified (09/07/2019) //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _MESHING_GEOMETRICAL_METHODS_H_
#define _MESHING_GEOMETRICAL_METHODS_H_


#include "MeshingCommon.h"
#include "MeshingRandomSampler.h"

class MeshingGeometricalMethods
{


public:

	MeshingGeometricalMethods();

    ~MeshingGeometricalMethods();

	double get_point_angle(double x, double y);

	double distance(size_t num_dim, double* x, double* y);

	double cos_angle(size_t num_dim, double* x, double* y, double* z);

	void get_power_vertex(size_t num_dim, double* co, double ro, double* c1, double r1, double* pv);

	bool get_power_vertex(size_t num_dim, size_t num_points, double** centers, double* radii, double* pv);

	int get_normal_component(size_t num_dim, size_t num_basis, double** basis, double* vect, double &norm);

	double dot_product(size_t num_dim, double* v1, double* v2);

	bool normalize_vector(size_t num_dim, double* vec);

	bool overlapping_spheres(size_t num_dim, double* sphere_i, double* sphere_j);

	double get_3d_triangle_area(double** corners);

	bool get_3d_triangle_normal(double** corners, double* normal);

	int project_to_3d_triangle(double* p, double** corners, double* q, double &proj_dist);

	int project_to_3d_line_segment(double* p, double* xo, double* xn, double* q, double& proj_dist);

	int project_to_3d_line(double* p, double* xo, double* edir, double* q, double &proj_dist);

private:
	MeshingRandomSampler _rsampler;
};

#endif

