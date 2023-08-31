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
//  MeshingVoroSpokes.h                                           Last modified (06/06/2019) //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _VORO_SPOKES_H_
#define _VORO_SPOKES_H_

#ifndef DBL_MAX
#define DBL_MAX          1.7976931348623158e+308 // max value
#endif


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
#include "MeshingRandomSampler.h"
#include "MeshingSmartTree.h"
#include "MeshingLinearSolver.h" // JHS

class MeshingVoroSpokes
{
	enum basis{monomials, Chebyshev};

public:

	MeshingVoroSpokes();

	~MeshingVoroSpokes();

	int init_vorospokes(size_t num_threads, size_t num_dim, size_t num_functions, size_t desired_order, double* xmin, double* xmax, bool adaptive_sampling_mode);

	int init_vorospokes(size_t num_threads, size_t num_dim, size_t budget, size_t num_functions, size_t desired_order, double* xmin, double* xmax, bool adaptive_sampling_mode);

	int clear_memory();

	int plot_function_slice(std::string file_name, double* po, double slice_width, size_t dim_1, size_t dim_2, size_t function_index, size_t num_contours, double f_min, double f_max, bool use_vwn, bool plot_points);

        int add_samples(size_t num_samples, double** samples_x, double** samples_f, size_t num_spokes);

	int set_sampling_function(size_t function_index, size_t num_spokes);

	int suggest_new_sample(double* xsample, size_t num_spokes);

	int add_new_seed(double* xsample, double* fsample, size_t num_spokes);

	// int get_closest_seed(double* x, size_t &iclosest, double &hclosest);

	double* get_seed(size_t seed_index) { return _seeds->get_tree_point(seed_index); };

	int evaluate_vps(double* x, double* fx);

	int evaluate_vps(double* x, double* fx, double &vps_err_est);

	int evaluate_vwn(double* x, size_t num_spokes, size_t num_layers, double* fx, double &vwn_err_est);

	int plot_function(std::string file_name, size_t function_index, size_t num_contours);

	size_t get_num_basis(size_t num_dim, size_t order);

	double evaluate_basis(size_t ibasis, size_t num_dim, size_t order, double* x, size_t &basis_order);

	int CG_LS(size_t num_rows, size_t num_columns, double** A, double* b, size_t nit_max, double eps, double* x);

	int build_vps_model(size_t icell, size_t num_spokes);
private:

	int init_random_samplers(size_t num_threads);

	int expand_vorospokes_containers();

	double evaluate_basis(double x, size_t ibasis);

	double evaluate_basis_chebyshev(double x, double fm, double fmm, size_t ibasis, size_t& jbasis);

	int evaluate_vps(size_t icell, double* x, double* f_x, double &vps_err_est);

	int evaluate_vps(size_t icell, double* x, double* f_x);

	int build_vps_model(size_t icell, size_t num_basis, double** c);

	int get_vps_neighbors(double* x, size_t num_spokes, size_t &num_neighbors, size_t*& neighbors, size_t*& neighbor_hits);

	double spherical_integration_constant(size_t num_dim);

	int integrate_vps(size_t icell, size_t num_spokes, size_t function_index, double &integration_val, double* xsample);

	int integrate_vps(size_t icell, size_t num_spokes, size_t function_index, double &integration_val, double* xsample, size_t thread_id);



private:
	size_t      _num_dim;
	size_t      _num_seeds;
	size_t      _budget;
	size_t      _num_functions;
	size_t      _desired_order;
	size_t*     _desired_orders;
	size_t      _sampling_function_index;

	double* _xmin;
	double* _xmax;

	MeshingSmartTree * _seeds;
	double**    _f;
	size_t**    _cell_order;
	double*     _cell_weight;
	double*     _cell_cdf;
	double**    _cell_samples;
	double***   _cell_coef;

	bool        _use_total_order;
	basis       _vps_basis;

	MeshingRandomSampler _rsampler;
	std::vector<MeshingRandomSampler*> _rsamplers;
	MeshingLinearSolver	_lin;


	};

#endif

