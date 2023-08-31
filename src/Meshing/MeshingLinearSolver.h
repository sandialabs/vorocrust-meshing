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
//  MeshingLinearSolver.h                                         Last modified (09/05/2019) //
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _LINEAR_SOLVER_H_
#define _LINEAR_SOLVER_H_

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

class MeshingLinearSolver
{

public:

	MeshingLinearSolver();

	~MeshingLinearSolver();

	int LS_QR_Solver(size_t nrow, size_t ncol, double** A, double* b, double* x);

	int householder(size_t nrow, size_t ncol, double** A, double** Q, double** R);

	bool Cholesky_factorization(size_t num_rows, size_t** AJ, double** AV);

	int Cholesky_factorization_rank_one_update(size_t num_rows, size_t** CJ, double** CV, size_t* &xj, double* &xv);

	int Cholesky_factorization_remove_row(size_t &num_rows, size_t** &CJ, double** &CV, size_t row_index);

	bool Cholesky_factorization_add_row(size_t& num_rows, size_t**& CJ, double**& CV, size_t* aj, double* av);

	int Cholesky_Solver(size_t num_rows, size_t** CJ, double** CV, double* b, double* x);


private:

	int set_sparse_vector_entry(size_t row, double val, size_t*& I, double*& V);

	int get_sparse_vector_entry(size_t row, double& val, size_t* I, double* V);

	size_t find_entry_position_binary(size_t entry, size_t* I, size_t left, size_t right);

	double get_dot_product(size_t* I1, double* V1, size_t* I2, double* V2, size_t min_index, size_t max_ind);

	double get_dot_product(size_t* I1, double* V1, double* V2, size_t min_index, size_t max_ind);
	
};

#endif

