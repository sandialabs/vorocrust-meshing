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
// MeshingVoroCrustSampler.h                                      Last modified (01/15/2019) //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _VOROCRUST_SAMPLER_H_
#define _VOROCRUST_SAMPLER_H_

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
#include <iomanip>

#include <time.h>

#include "MeshingCommon.h"
#include "MeshingPolyMesh.h"
#include "MeshingSmartTree.h"
#include "MeshingRandomSampler.h"
#include "MeshingGeometricalMethods.h"
#include "MeshingMemoryHandler.h"
#include "MeshingSpheresMethods.h"
#include "MeshingVoroSpokes.h" // What is VoroSpokes?

class MeshingVoroCrustSampler
{

public:

	MeshingVoroCrustSampler(size_t num_threads);

    ~MeshingVoroCrustSampler();

	int set_input_sizing_function(MeshingSmartTree* sizing) { _input_sizing_function = sizing; return 0; };

	int dont_detect_features_generate_point_clouds(size_t num_points, double** points, size_t num_faces, size_t** faces,
		                                           MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud,
		                                           size_t& num_sharp_edges, size_t**& sharp_edges, size_t& num_sharp_corners, size_t*& sharp_corners);

	int detect_features_generate_point_clouds_clean(size_t num_points, double** points, size_t num_faces, size_t** faces,
		                                            double smooth_angle_threshold,
		                                            MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud,
		                                            size_t& num_sharp_edges, size_t**& sharp_edges, size_t& num_sharp_corners, size_t*& sharp_corners);


	int detect_features_generate_point_clouds(size_t num_points, double** points, size_t num_faces, size_t** faces,
		                                      double smooth_angle_threshold, double rmin,
		                                      MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud,
		                                      size_t &num_sharp_edges, size_t** &sharp_edges, size_t &num_sharp_corners, size_t* &sharp_corners);


	int generate_corner_spheres(size_t num_points, double** points, size_t num_sharp_corners, size_t* sharp_corners,
		                        MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud, double feature_TOL,
		                        MeshingSmartTree* corner_spheres,
		                        double smooth_angle_threshold, double rmax, double Lip);

	int generate_edge_spheres(size_t num_points, double** points, size_t num_edges, size_t** edges,
		                      MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud, double feature_TOL,
		                      MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres,
		                      double smooth_angle_threshold, double rmax, double Lip, double coverage_radius_ratio);


	int generate_surface_spheres(size_t num_points, double** points, size_t num_faces, size_t** faces,
		                         MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud, double rmin,
		                         MeshingSmartTree* &surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres,
		                         MeshingSmartTree* &ext_surf_seeds, MeshingSmartTree* &int_surf_seeds,
		                         double smooth_angle_threshold, double rmax, double Lip, double coverage_radius_ratio);

	int color_surface_seeds(MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres, MeshingSmartTree* upper_seeds, MeshingSmartTree* lower_seeds,
		                    MeshingSmartTree* spheres, MeshingSmartTree* seeds, size_t &num_subregions);
	
private:

	int clear_memory();

	// VC WILD
	int initiate_active_pool_edges(size_t num_points, double** points, size_t num_edges, size_t** edges);

	// VC WILD
	int refine_active_pool_edges(size_t num_points, double** points, size_t num_edges, size_t** edges,
		                         MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres, double Lip, double coverage_radius_ratio);

	int sample_active_edges_uniformly(size_t num_points, double** points, size_t num_edges, size_t** edges, double* dart, size_t &dart_edge);

	int get_active_edge_corners(size_t num_points, double** points, size_t num_edges, size_t** edges, size_t active_edge_index, double** edge_corners);

	int refine_active_edge(size_t active_edge_index, double** edge_corners,
		                   size_t &active_pool_size, size_t &active_pool_capacity,
		                   size_t* &active_pool_parent_face, bool** &active_pool_children, double* &active_pool_cdf);

	// VC WILD
	int initiate_active_pool(size_t num_points, double** points, size_t num_faces, size_t** faces, size_t active_patch, size_t &num_patches);

	int sample_active_pool(size_t num_samples, size_t num_points, double** points, size_t num_faces, size_t** faces,
		                   MeshingSmartTree* surface_point_cloud, MeshingSmartTree* edge_point_cloud, double R_MIN,
		                   MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres,
		                   double smooth_angle_threshold, double rmax, double Lip, double coverage_radius_ratio,
	                       double** new_sphere, size_t* new_sphere_face, bool* covered_sample);

	int refine_active_pool(size_t num_points, double** points, size_t num_faces, size_t** faces,
		                   MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres, double Lip, double coverage_radius_ratio);

	int sample_active_faces_uniformly(size_t thread_id, size_t num_points, double** points, size_t num_faces, size_t** faces, double* dart, size_t& dart_face);

	int get_active_face_corners(size_t num_points, double** points, size_t num_faces, size_t** faces, size_t active_face_index, double** face_corners);

	int refine_active_face(size_t active_face_index, double** face_corners,
		                   size_t &active_pool_size, size_t &active_pool_capacity,
		                   size_t* &active_pool_parent_face, bool** &active_pool_children, double* &active_pool_cdf);

	int delete_active_pool();

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int generate_surface_seeds(size_t num_points, double** points, size_t num_faces, size_t** faces,
		                       MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres,
		                       double Lip, bool &shrunk_s,
		                       MeshingSmartTree* upper_seeds, MeshingSmartTree* lower_seeds);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	bool validate_surface_mesh(MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres, size_t num_faces, size_t** faces);

private:

	double     _ymin;

	double     _r_min_edges;
	double     _r_min_surface;

	size_t     _num_faces;
	double**   _face_normal;


	size_t	   _active_pool_ref_level;
	size_t     _active_pool_size;
	size_t*    _active_pool_parent_face;
	bool**     _active_pool_children;
	double*    _active_pool_cdf;

	double _pool_sampling_time;
	double _pool_refinement_time;

	double _sizing_time;

	MeshingRandomSampler				_rsampler;

	size_t                      _num_threads;
	std::vector<MeshingRandomSampler*> _rsamplers;


	MeshingGeometricalMethods		_geom;
	MeshingMemoryHandler           _memo;
	MeshingSpheresMethods          _smethods;
	MeshingPolyMesh                _mesh;

	MeshingVoroSpokes              _surface_sizing_function;

	MeshingSmartTree*              _input_sizing_function;
};

#endif

