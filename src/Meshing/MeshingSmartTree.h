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
//  MeshingSmartTree.h                                            Last modified (07/12/2022) //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _MESHING_SMART_TREE_H_
#define _MESHING_SMART_TREE_H_


#include "MeshingCommon.h"
#include "MeshingLinearSolver.h"
#include "MeshingRandomSampler.h"
#include "MeshingGeometricalMethods.h"
#include "MeshingMemoryHandler.h"

class MeshingSmartTree
{

public:

	MeshingSmartTree();

	MeshingSmartTree(size_t num_dim);

	MeshingSmartTree(size_t num_dim, size_t num_points, double** points);

	~MeshingSmartTree();


	int re_enumerate_points_for_better_memory_access();

	int init_random_samplers(size_t num_threads);

	int disable_auto_balance() { _auto_balance = false; return 0;}

	int enable_auto_balance() { _auto_balance = true; return 0; }

	int save_tree_csv(std::string file_name, size_t num_dim);

	int get_tokens(std::string line, char separator, std::vector<std::string>& tokens);

	int load_tree_csv(std::string file_name, size_t num_dim);

	int reset_graph();

	int graph_get_neighbors(size_t ipoint, size_t*& neighbors) { neighbors = _graph[ipoint]; return 0;};

	// connect non-directional edge
	int graph_connect_nodes(size_t ipoint, size_t jpoint);

	// connect directional_edge
	int graph_connect(size_t ipoint, size_t jpoint);

	// disconnect non-directional edge
	int graph_disconnect_nodes(size_t ipoint, size_t jpoint);

	// disconnect directional edge
	int graph_disconnect(size_t ipoint, size_t jpoint);

	bool graph_connected(size_t ipoint, size_t jpoint);

	size_t get_num_dimensions() { return _num_dim; };

	size_t get_num_tree_points() { return _num_points; };

	int get_bounding_box(size_t &num_dim, double* &xmin, double* &xmax);

	void set_tree_points(size_t num_dim, size_t num_points, double** points);

	int add_tree_point(size_t num_dim, double* x, double* normal, size_t* attrib);

	int get_tree_point(size_t point_index, double* x);

	int get_tree_point(size_t point_index, size_t num_dim, double* x);

	bool get_tree_point_attrib(size_t point_index, size_t attrib_index, size_t &point_attrib);

	double get_tree_point_attrib(size_t point_index, size_t attrib_index);

	double* get_tree_point(size_t point_index) {
		return _points[point_index];
	};

	double* get_tree_point_normal(size_t point_index) {
		return _points_normal[point_index];
	};

	size_t* get_tree_point_attrib(size_t point_index) {
		return _points_attrib[point_index];
	};

	int set_tree_point_attrib(size_t point_index, size_t attrib_index, size_t attrib);

	int set_tree_point_attrib(size_t point_index, size_t attrib_index, double attrib);

	void set_marked_tree_points(bool* marked) { _marked_only = true; _marked = marked; };

	void deactivate_marked_tree_points() { _marked_only = false;};

	void clear_memory();

	int kd_tree_build_balanced();

	int get_closest_tree_point(size_t tree_point_index, size_t &closest_tree_point, double &closest_distance);

	int get_closest_tree_point(double* x, size_t &closest_tree_point, double &closest_distance);

	int get_closest_tree_point(double* x, size_t max_index, size_t &closest_tree_point, double &closest_distance);

	int get_closest_tree_point(double* x, size_t &closest_tree_point, double &closest_distance, size_t &num_nodes_visited);

	int get_closest_tree_point(double* x, size_t num_exculded_tree_points, size_t* exculded_tree_points, size_t &closest_tree_point, double &closest_distance);

	int get_closest_non_smooth_tree_point(bool directional_normal, double* x, size_t num_xnormals, double** xnormals, double smoothness_threshold, size_t &closest_tree_point, double &closest_distance);

	// used for estimating surface sizing function in VC_Wild
	int get_closest_non_smooth_tree_point(double* x, double feature_TOL, double smoothness_threshold, size_t &closest_tree_point, double &closest_distance);

	// used for estimating surface sizing function in VC_Wild
	int get_closest_non_smooth_tree_point(double* x, double* x_normal, double smoothness_threshold, size_t& closest_tree_point, double& closest_distance);

	// used for detecting sharp corners/edges in VC_Wild
	int get_closest_non_smooth_tree_point(double* x, size_t num_normals, double** normals, double smoothness_threshold, size_t &closest_tree_point, double &closest_distance);

	// used for estimating curves sizing function in VC_Wild
	int get_closest_non_smooth_tree_edge_point(double* x, double feature_TOL, double smoothness_threshold, size_t &closest_tree_point, double &closest_distance);

	// used for detecting sharp corners due to edges in VC_Wild
	int get_closest_non_smooth_tree_edge_point(double* x, double* tangent, double smoothness_threshold, size_t &closest_tree_point, double &closest_distance);



	int get_tree_points_in_sphere(double* x, double r, size_t &num_points_in_sphere, size_t* &points_in_sphere, size_t &capacity);

	int get_tree_points_in_sphere(double* x, double r, size_t max_index, size_t &num_points_in_sphere, size_t* &points_in_sphere, size_t &capacity);

	int get_tree_points_in_spheres(size_t num_in_spheres, double** in_spheres, double* r_in,
		                           size_t num_out_spheres, double** out_spheres, double* r_out,
		                           size_t &num_points_in_sphere, size_t* &points_in_sphere, size_t &capacity);

	int get_tree_point_in_spheres(size_t num_in_spheres, double** in_spheres, double* r_in,
		                          size_t num_out_spheres, double** out_spheres, double* r_out, size_t &tree_point_index);

	int trim_spoke(double* cell_seed, size_t num_involved_seeds, size_t* involved_seeds, double* spoke_origin, double* spoke_dir, double& spoke_length, size_t& neighbor_seed);

	int get_Voronoi_neighbors(size_t num_spokes, size_t num_layers, size_t num_dim, double* x,
		                      size_t &num_neighbors, size_t* &neighbors, size_t* &neighbor_hits, size_t* &neighbor_level);

	int get_Voronoi_neighbors(size_t num_spokes, size_t num_layers, size_t num_dim, double* x,
		                      double** spoke_points, size_t* spoke_neighbors,
		                      size_t &num_neighbors, size_t* &neighbors, size_t* &neighbor_hits);

	int get_Voronoi_neighbors(size_t num_desired_neighbors, size_t num_dim, double* x,
		                      size_t& num_neighbors, size_t*& neighbors);

	int get_Voronoi_neighbors(size_t num_spokes, size_t num_layers, size_t num_dim, double* x,
		                      double** spoke_points, size_t* spoke_neighbors,
		                      size_t &num_neighbors, size_t* &neighbors, size_t* &neighbor_hits, size_t* &neighbor_level);

	int get_Voronoi_neighbors(size_t thread_id, size_t num_spokes, size_t num_layers, size_t num_dim, double* x,
		                      double** spoke_points, size_t* spoke_neighbors,
		                      size_t& num_neighbors, size_t*& neighbors, size_t*& neighbor_hits);

	int get_Voronoi_neighbors(size_t thread_id, size_t num_spokes, size_t num_layers, size_t num_dim, double* x,
		                      double** spoke_points, size_t* spoke_neighbors,
		                      size_t& num_neighbors, size_t*& neighbors, size_t*& neighbor_hits, size_t*& neighbor_level);

	int get_3d_Voronoi_cell(size_t seed_index, size_t &num_cell_neighbors, size_t* &cell_neighbors, double &cell_volume,
		                    size_t* &num_faces_corners, double*** &faces_corners, double* &faces_area);

	int get_3d_Voronoi_face(size_t iseed, size_t jseed, double* xface,
		                    size_t &num_face_corners, double** &face_corners,
		                    double &face_area, double &face_height);

	int get_2d_Voronoi_cell(size_t iseed, double* xmin, double* xmax, size_t& num_corners, double** &corners, size_t* &neighbors);

	int get_2d_Voronoi_cell(double* x, double* xmin, double* xmax, size_t &num_corners, double** &corners, size_t* &neighbors);

	bool get_Voronoi_vertex(size_t num_seeds, size_t* seeds, double* vertex);

	int plot_delaunay_graph(std::string outFile, size_t num_spokes, size_t num_layers);

	int plot_points(std::string outFile, bool verbose);

	int save_Voronoi_face_ply(std::string file_name, size_t num_face_corners, double** face_corners);

	int save_Voronoi_faces_ply(std::string file_name, size_t num_faces, size_t* num_face_corners, double*** face_corners);

	int save_Voronoi_faces_ply(std::string file_name, size_t num_faces, size_t* num_face_corners, double*** face_corners, double xmax, double ymax, double zmax);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////// kd-tree Private Methods       /////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

private:

	int init_global_variables();

	double string_to_double(const std::string &s)
	{
		std::istringstream i(s); double x;
		if (!(i >> x)) return 0.0; return x;
	};

	size_t string_to_size_t(const std::string &s)
	{
		std::istringstream i(s); size_t x;
		if (!(i >> x)) return 0; return x;
	};

	int get_normal_component(size_t num_dim, size_t num_basis, double** basis, double* vect, double &norm);

	int kd_tree_balance_quicksort(size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted);

	int kd_tree_quicksort_adjust_target_position(size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted);

	int kd_tree_add_point(size_t point_index);

	int kd_tree_get_nodes_order(size_t d_index, size_t node_index,                                // indices to traverse the kd-tree
		                        size_t &num_traversed, size_t* ordered_indices
	                            );


	int kd_tree_get_seeds_in_sphere(double* x,                                                    // Sphere center
		                            double r,                                                     // Sphere radius
		                            size_t d_index, size_t node_index,                            // indices to traverse the kd-tree
		                            size_t &num_points_in_sphere, size_t* &points_in_sphere,      // number of points in sphere and their indices
		                            size_t &capacity                                              // Size of points in sphere array
	                               );

	int kd_tree_get_seeds_in_sphere(double* x,                                                    // Sphere center
		                            double r,                                                     // neighborhood radius
		                            size_t max_index,                                             // maximum index to be considered
		                            size_t d_index, size_t node_index,                            // indices to traverse the kd-tree
		                            size_t &num_points_in_sphere, size_t* &points_in_sphere,      // number of points in sphere and their indices
		                            size_t &capacity                                              // Size of points in sphere array
	                               );

	int kd_tree_get_seed_in_spheres(size_t num_in_spheres, double** in_spheres, double* r_in,
		                            size_t num_out_spheres, double** out_spheres, double* r_out,
		                            size_t d_index, size_t node_index, size_t &seed_index);

	int kd_tree_get_seeds_in_spheres(size_t num_in_spheres, double** in_spheres, double* r_in,
		                             size_t num_out_spheres, double** out_spheres, double* r_out,
		                             size_t d_index, size_t node_index,
		                             size_t &num_points_in_sphere, size_t* &points_in_sphere, size_t &capacity);

	int kd_tree_get_closest_seed(size_t seed_index,                                  // seed index
		                         size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
		                         size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
	                            );

	int kd_tree_get_closest_seed(double* x,                                          // a point in space
		                         size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
		                         size_t &closest_seed, double &closest_distance,      // index of closest seed and distance from it
		                         size_t &num_nodes_visited
	                            );

	int kd_tree_get_closest_seed(double* x,                                          // point location
		                         size_t max_index,                                   // maximum index to be considered
		                         size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
		                         size_t &closest_seed, double &closest_distance,      // index of closest seed and distance from it
		                         size_t &num_nodes_visited
	                            );

	int kd_tree_get_closest_seed(double* x,                                          // a point in space
		                         size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
		                         size_t num_exculded_seeds, size_t* exculded_seeds,  // Number of excluded seeds and their indices
		                         size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
	                            );

	int kd_tree_get_closest_non_smooth_seed(bool directional_normal,                            // use directional point normal
		                                    double* x,                                          // point location
		                                    size_t num_normals, double** x_normals,             // point normals
		                                    double smoothness_threshold,                        // smoothness threshold
		                                    size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
		                                    size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
	                                       );

	int kd_tree_get_closest_non_smooth_seed(double* x, double feature_TOL,                  // point location
		                                    double smoothness_threshold,                    // smoothness threshold
		                                    size_t d_index, size_t node_index,              // indices to traverse the kd-tree
		                                    size_t &closest_seed, double &closest_distance  // index of closest seed and distance from it
	                                      );

	int kd_tree_get_closest_non_smooth_seed(double* x, double* x_normal,                    // point location and normal
		                                    double smoothness_threshold,                    // smoothness threshold
		                                    size_t d_index, size_t node_index,              // indices to traverse the kd-tree
		                                    size_t& closest_seed, double& closest_distance  // index of closest seed and distance from it
	                                       );

	int kd_tree_get_closest_non_smooth_seed(double* x, size_t num_normals, double** normals,                      // point location and normal
		                                    double smoothness_threshold,                                          // smoothness threshold
		                                    size_t d_index, size_t node_index,                                    // indices to traverse the kd-tree
		                                    size_t &closest_seed, double &closest_distance                        // index of closest seed and distance from it
	                                        );

	int kd_tree_get_closest_non_smooth_edge_seed(double* x, double feature_TOL,                   // point location and feature tolernace
		                                         double smoothness_threshold,                    // smoothness threshold
		                                         size_t d_index, size_t node_index,              // indices to traverse the kd-tree
		                                         size_t &closest_seed, double &closest_distance  // index of closest seed and distance from it
	                                            );

	int kd_tree_get_closest_non_smooth_edge_seed(double* x, double* tangent,                      // point location and tangent
		                                         double smoothness_threshold,                    // smoothness threshold
		                                         size_t d_index, size_t node_index,              // indices to traverse the kd-tree
		                                         size_t &closest_seed, double &closest_distance  // index of closest seed and distance from it
	                                            );

	int kd_tree_trim_spoke(double* cell_seed, // could be a new seed(no in _points)
		                   size_t num_involved_seeds, size_t* involved_seeds, double* spoke_origin, double* spoke_dir, double origin_dst,
		                   size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
		                   double &spoke_length, size_t &neighbor_seed         // spoke length and last trimming neighbor
	                      );

	int get_closest_seed_brute(double* x, size_t num_significant_neighbors, size_t* significant_neighbors, size_t num_exculded_seeds, size_t* exculded_seeds, size_t &closest_seed);

	int add_entry(size_t entry, size_t &num_entries, size_t* &I, size_t &capacity);

	int add_entry(double* entry, size_t &num_entries, double** &I, size_t &capacity);

	bool find_brute(size_t entry, size_t* I, size_t num_entries);

	int quicksort(double* x, double* y1, double* y2, size_t* I, size_t left, size_t right);

	bool find_face_vertex(double* si, double* sj, double* qm, double* qp, double** e_basis, double &tm, double &tp, double &vx, double &vy);

	////////////////////////////////////////////////////////////////////////////////
	/////////// kd-tree  Variables /////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	size_t            _num_dim;                    // number of dimensions
	size_t            _capacity;
	size_t            _num_points;                 // number of points
	double**          _points;                     // points
	double**          _points_normal;              // points_normal

	double*           _xmin;
	double*           _xmax;

	size_t**          _points_attrib;              // a list of integers asscoaited with each point

	size_t            _graph_cap;           // size of graph list
	size_t**          _graph;               // a list of integers asscoaited with each point

	size_t            _tree_origin;                // index of tree root
	size_t            _tree_max_height;            // Max height of a tree branch
	size_t*           _tree_left;                  // left pointer of a tree node
	size_t*           _tree_right;                 // right pointer of a tree node

	bool              _auto_balance;

	bool             _marked_only;                // default is false
	bool*            _marked;

	MeshingRandomSampler               _rsampler;
	std::vector<MeshingRandomSampler*> _rsamplers;
	MeshingLinearSolver                _lin_solver;
	MeshingGeometricalMethods          _geom;
	MeshingMemoryHandler               _memo;

};

#endif

