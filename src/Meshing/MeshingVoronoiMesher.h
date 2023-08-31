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
//  MeshingVoronoiMesher.h                                        Last modified (11/02/2020) //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MeshingCommon.h" // JHS 08/03/2022 changed from Common.h to MeshingCommon.h
#include "MeshingSmartTree.h"
#include "MeshingTimer.h"
#ifndef _VORONOI_MESHER_H_
#define _VORONOI_MESHER_H_

class MeshingVoronoiMesher
{

public:

	MeshingVoronoiMesher();

	~MeshingVoronoiMesher() { };

	int ensure_sharp_edge_spheres_are_Delaunay(size_t& num_seeds, double*& seeds, double*& seeds_sizing, size_t*& seeds_region_id,
		                                       size_t num_surface_spheres, double* surface_spheres, double* surface_spheres_sizing,
		                                       size_t num_edge_spheres, double* edge_spheres, double* edge_spheres_sizing,
		                                       size_t num_corner_spheres, double* corner_spheres, double* corner_spheres_sizing,
		                                       size_t num_sizing_points, double* sizing_points, double* sizing_value,
		                                       int num_threads, double Lip, double rmax);

	int impose_interior_seeds(size_t& num_seeds, double*& seeds, double*& seeds_sizing, size_t*& seeds_region_id,
		                      size_t num_surface_spheres, double* surface_spheres, double* surface_spheres_sizing,
		                      size_t num_edge_spheres, double* edge_spheres, double* edge_spheres_sizing,
		                      size_t num_corner_spheres, double* corner_spheres, double* corner_spheres_sizing,
		                      size_t num_sizing_points, double* sizing_points, double* sizing_value,
		                      size_t num_imposed_seeds, double* imposed_seeds,
		                      int num_threads, double Lip, double rmax);

	int generate_interior_seeds(size_t& num_seeds, double*& seeds, double*& seeds_sizing, size_t*& seeds_region_id,
		                        size_t num_surface_spheres, double* surface_spheres, double* surface_spheres_sizing,
		                        size_t num_edge_spheres, double* edge_spheres, double* edge_spheres_sizing,
		                        size_t num_corner_spheres, double* corner_spheres, double* corner_spheres_sizing,
		                        size_t num_sz_points, double* sz_points, double* sz_value,
		                        int num_threads, double Lip, double rmax);

	int generate_3d_voronoi_mesh(int num_threads, size_t num_seeds, double* seeds, size_t* seed_region_id, double* seed_sizing,
		                         size_t& num_vertices, double*& vertices, size_t& num_faces, size_t** &faces);

	int get_3d_Voronoi_cell(size_t num_seeds, double* seeds, size_t tree_origin, size_t* tree_right, size_t* tree_left,
		                    size_t seed_index, double r, bool* ghost_seed,
		                    size_t& num_cell_corners, double*& cell_corners, size_t*& cell_corner_neighbors,
		                    size_t*& cell_corner_seeds, size_t*& cell_corner_indices);

	int save_vcg_file_pflotran(std::string outFile, size_t num_seeds, double* seeds, size_t* seed_region_id,
		                       size_t num_vertices, double* vertices, size_t num_faces, size_t** faces,
		                       size_t* cell_num_faces, size_t** cell_faces);

    int save_exodus_file(std::string file_name, size_t num_seeds, double* seeds, size_t* seed_region_id,
                         size_t num_vertices, double* vertices, size_t num_faces, size_t** faces,
                         size_t* cell_num_faces, size_t** cell_faces);
    
	int save_spheres_csv(std::string file_name, size_t num_dim, size_t num_spheres, 
		                 double* spheres, double* spheres_sizing, size_t* spheres_region_id);

	int load_spheres_csv(std::string file_name, size_t num_dim, size_t& num_spheres, double*& spheres);

	int save_Voronoi_tessellation(std::string file_name, size_t num_vertices, double* vertices, size_t num_faces, size_t** faces);

	int save_Voronoi_tessellation(std::string file_name, size_t num_vertices, double* vertices, size_t num_faces, size_t** faces, double* xo, double* edir);

private:
	int construct_initial_polytope(double* x, double r, 
		                           size_t& num_corners, size_t& corners_cap, double*& corners, size_t*& corners_neighbors, 
		                           size_t*& corners_seeds);

	int trim_vertex(size_t trimming_seed_index, double* xo, double* n, size_t corner_index, size_t old_corner,
		            size_t& num_corners, size_t& corners_cap, double*& corners, size_t*& corners_neighbors, size_t*& corners_seeds);

	
	int build_balanced_kd_tree(size_t num_points, size_t num_dim, double* points, size_t& tree_origin, size_t* tree_right, size_t* tree_left);

	int re_enumerate_points_for_better_memory_access(size_t num_points, size_t num_dim, double* points,
		                                             size_t& tree_origin, size_t* tree_right, size_t* tree_left,
		                                             double* points_sorted, size_t* point_old_index, size_t* point_new_index);

	int get_closest_tree_point(size_t num_points, size_t num_dim, double* points, size_t tree_origin, size_t* tree_right, size_t* tree_left,
		                       double* x, size_t& closest_tree_point, double& closest_distance);

	int get_closest_tree_point(size_t num_points, size_t num_dim, double* points, size_t tree_origin, size_t* tree_right, size_t* tree_left,
		                       size_t tree_point_index, double* e_dir, size_t& closest_tree_point, double& closest_distance);


	int get_tree_points_in_sphere(size_t num_points, size_t num_dim, double* points, size_t tree_origin, size_t* tree_right, size_t* tree_left,
		                          double* x, double r, size_t& num_points_in_sphere, size_t*& points_in_sphere);

	int kd_tree_balance_quicksort(size_t num_points, size_t num_dim, double* points,
		                          size_t &tree_origin, size_t* tree_right, size_t* tree_left, size_t& tree_height,
		                          size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted);

	int kd_tree_quicksort_adjust_target_position(size_t num_points, size_t num_dim, double* points, 
		                                         size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted);


	int kd_tree_add_point(size_t num_points, size_t num_dim, double* points,
		                  size_t tree_origin, size_t* tree_right, size_t* tree_left, size_t& tree_height,
		                  size_t seed_index);

	int kd_tree_get_closest_seed(size_t num_points, size_t num_dim, double* points,
		                         size_t tree_origin, size_t* tree_right, size_t* tree_left,
		                         double* x, size_t d_index, size_t node_index,
		                         size_t& closest_seed, double& closest_distance,
		                         size_t& num_nodes_visited);

	int kd_tree_get_closest_seed(size_t num_points, size_t num_dim, double* points,
		                         size_t tree_origin, size_t* tree_right, size_t* tree_left,
		                         size_t tree_point_index, double* e_dir,
		                         size_t d_index, size_t node_index,
		                         size_t& closest_seed, double& closest_distance,
		                         size_t& num_nodes_visited);

	int kd_tree_get_seeds_in_sphere(size_t num_points, size_t num_dim, double* points,
		                            size_t tree_origin, size_t* tree_right, size_t* tree_left,
		                            double* x, double r, size_t d_index, size_t node_index,
		                            size_t& num_points_in_sphere, size_t*& points_in_sphere, size_t& capacity);

	int kd_tree_get_nodes_order(size_t* tree_right, size_t* tree_left, 
		                        size_t d_index, size_t node_index,
		                        size_t& num_traversed, size_t* ordered_indices);

	
	int impose_lipschitz_continuity(size_t num_spheres, size_t num_dim, double* spheres, double* spheres_sizing, 
		                            int num_threads, double Lip,
		                            size_t tree_origin, size_t* tree_right, size_t* tree_left);

	int get_tokens(std::string line, char separator, std::vector<std::string>& tokens);

	MeshingMemoryHandler	_memo;

    
   

};

#endif

