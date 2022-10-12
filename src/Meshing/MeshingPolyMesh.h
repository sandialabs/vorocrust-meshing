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
// MeshingPolyMesh.h                                              Last modified (01/12/2018) //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _POLY_MESH_H_
#define _POLY_MESH_H_


#include "MeshingCommon.h"
#include "MeshingRandomSampler.h"
#include "MeshingSpheresMethods.h"
#include "MeshingGeometricalMethods.h"

#include "MeshingMemoryHandler.h"
#include "MeshingSmartTree.h"

class MeshingPolyMesh
{

public:

	MeshingPolyMesh();

	MeshingPolyMesh(size_t num_points, double** points, size_t num_faces, size_t** faces, double smooth_angle_threshold);

    ~MeshingPolyMesh();

	int clear_memory();

	int read_input_obj_file(std::string filename, double smooth_angle_threshold);

	int read_input_obj_file(std::string filename, size_t &num_points, double** &points, size_t &num_faces, size_t** &faces);

	int save_mesh_obj(std::string file_name, size_t num_points, double** points, size_t num_faces, size_t** faces);

	int weld_nearby_points(size_t &num_points, double** &points, size_t &num_faces, size_t** &faces, double hmin);

	bool project_point_to_neighbors_convex_hull(double* x, size_t num_neighbors, double** neighbors, double* proj);

	int process_skinny_triangles(size_t num_points, double** points, size_t num_faces, size_t** faces);

	int extract_connected_smooth_layers();

	bool generate_surface_mesh(MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres,
		                       double Lip, bool &shrunk_s, bool &shrunk_e, bool &shrunk_c);

	bool point_covered(double* point, MeshingSmartTree* spheres, double Lip, double alpha_coverage, size_t si, size_t sj, size_t sk,
		               size_t &num_covering_spheres, size_t& cap_covering_spheres, size_t* &covering_spheres);

	bool valid_surface_mesh(MeshingSmartTree* surface_spheres, MeshingSmartTree* edge_spheres, MeshingSmartTree* corner_spheres, size_t num_faces, size_t** faces, bool &surface_sphere_shrunk, bool &edge_sphere_shrunk, bool &corner_sphere_shrunk);

	int process_model();

	bool water_tight_model();

	bool water_tight_model(size_t num_points, double** points, size_t num_faces, size_t** faces);

	int impose_consistent_orientation();

	int extract_model_sharp_features();

	int build_corners_point_cloud();

	int build_edge_point_cloud(size_t num_points);

	int build_surface_point_cloud(size_t num_points);

	int identify_face_smooth_neighbors();

	int identify_sharp_edge_smooth_neighbors();

	int improve_aspect_ratio(size_t num_refinements);

	int smooth_model_via_loop_subdivision(size_t num_refinements);

	size_t get_num_edge_faces(size_t point_i, size_t point_j);

	size_t get_num_edge_faces(size_t point_i, size_t point_j, size_t** faces, size_t** point_faces);

	int get_edge_faces(size_t point_i, size_t point_j, size_t* &edge_faces);

	int get_edge_faces(size_t point_i, size_t point_j, size_t** faces, size_t** point_faces, size_t* &edge_faces);

	bool smooth_neighbor_faces(size_t face_i, size_t face_j);

	bool smooth_neighbor_sharp_edges(size_t edge_i, size_t edge_j);

	bool sharp_corner(size_t point_index);

	bool sharp_edge(size_t pi, size_t pj);

	int get_opposite_corner(size_t corner, size_t face_index, size_t &opposite_corner, size_t &opposite_face);

	int build_point_faces(size_t num_points, size_t num_faces, size_t** faces, size_t** &point_faces);

	int build_face_normals(double** points, size_t num_faces, size_t** faces, double** &face_normal);

	int remove_isolated_points(size_t &num_points, double** &points, size_t num_faces, size_t** faces);

private:

    int reset_global_variables();

	double string_to_double(const std::string &s);

	int remove_isolated_points();

	int build_bounding_box();

	int build_point_sharp_edges();

	int build_point_faces();

	int build_face_normals();

	int adjust_cad_orientation();

	int write_spheres_csv(std::string file_name, size_t num_spheres, double** spheres);


public:

	double      _diag;
	double*     _xmin;
	double*     _xmax;

	size_t	    _num_points;
	double**    _points;

	size_t      _num_faces;
	size_t**    _faces;

	size_t**    _point_faces;
	double**    _face_normal;
	size_t**    _face_neighbors;

	size_t      _num_sharp_corners;
	size_t*     _sharp_corners;

	size_t      _num_sharp_edges;
	size_t**    _sharp_edges;
	size_t**    _point_sharp_edges;

	size_t**    _sharp_edge_neighbors;

	double**    _point_basis;

	double      _smooth_angle_threshold;

	MeshingGeometricalMethods      _geom;
	MeshingRandomSampler           _rsampler;
	MeshingSpheresMethods                 _smethods;
	MeshingMemoryHandler      _memo;

	MeshingSmartTree          _plc_corners_point_cloud;
	MeshingSmartTree          _plc_edge_point_cloud;
	MeshingSmartTree          _plc_surface_point_cloud;
};

#endif

