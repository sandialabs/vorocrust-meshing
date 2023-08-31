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
//  MeshingLinearSolver.cpp                                       Last modified (09/05/2019) //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MeshingLinearSolver.h"

#define SPARSE_ENTRY_TOL 1E-12

MeshingLinearSolver::MeshingLinearSolver()
{

}

MeshingLinearSolver::~MeshingLinearSolver()
{

}

int MeshingLinearSolver::LS_QR_Solver(size_t nrow, size_t ncol, double** A, double* b, double* x)
{
	#pragma region Least Square QR Solver:
	double** Q = new double* [nrow];
	for (size_t irow = 0; irow < nrow; irow++) Q[irow] = new double[nrow];

	double** R = new double* [nrow];
	for (size_t irow = 0; irow < nrow; irow++) R[irow] = new double[ncol];

	householder(nrow, ncol, A, Q, R);

	double* bmod = new double[ncol];
	for (size_t irow = 0; irow < ncol; irow++)
	{
		bmod[irow] = 0.0;
		for (size_t icol = 0; icol < nrow; icol++) bmod[irow] += Q[icol][irow] * b[icol];
	}

	// backward solve
	size_t jrow = ncol - 1;
	while (true)
	{
		if (fabs(R[jrow][jrow]) > 1E-10)
		{
			double sum(0.0);
			for (size_t jcol = ncol - 1; jcol > jrow; jcol--) sum += R[jrow][jcol] * x[jcol];
			x[jrow] = (bmod[jrow] - sum) / R[jrow][jrow];
		}
		else
		{
			// This column is in the Null Space of A
			x[jrow] = 0.0;
		}

		if (jrow == 0) break;
		jrow--;
	}

	for (size_t irow = 0; irow < nrow; irow++) delete[] Q[irow];
	delete[] Q;

	for (size_t irow = 0; irow < nrow; irow++) delete[] R[irow];
	delete[] R;

	delete[] bmod;
	return 0;
	#pragma endregion
}

int MeshingLinearSolver::householder(size_t nrow, size_t ncol, double** A, double** Q, double** R)
{
	#pragma region Householder QR Factorization:
	if (ncol > nrow)
	{
		vcm_cout << "Error (householder): Matrix has more columns than rows!" << std::endl;
		return 1;
	}

	for (size_t i = 0; i < nrow; i++)
	{
		for (size_t j = 0; j < ncol; j++) R[i][j] = A[i][j];
		for (size_t j = 0; j < nrow; j++) Q[i][j] = 0.0;
		Q[i][i] = 1.0;
	}

	double* x = new double[nrow];
	double* x_R = new double[ncol];
	double* x_Q = new double[nrow];
	for (size_t k = 0; k < ncol; k++)
	{
		double alpha(0.0);
		for (size_t i = k; i < nrow; i++)
		{
			x[i] = R[i][k];
			alpha += x[i] * x[i];
		}
		alpha = sqrt(alpha);
		if (x[k] > 0) alpha = -alpha;
		x[k] -= alpha;

		double norm(0.0);
		for (size_t i = k; i < nrow; i++) norm += x[i] * x[i];
		norm = sqrt(norm);
		for (size_t i = k; i < nrow; i++) x[i] /= norm;

		// update x_R
		for (size_t j = k; j < ncol; j++)
		{
			x_R[j] = 0.0;
			for (size_t i = k; i < nrow; i++)
			{
				x_R[j] += x[i] * R[i][j];
			}
		}

		// update x_Q
		for (size_t j = 0; j < nrow; j++)
		{
			x_Q[j] = 0.0;
			for (size_t i = k; i < nrow; i++)
			{
				x_Q[j] += x[i] * Q[i][j];
			}
		}

		// update R
		for (size_t i = k; i < nrow; i++)
		{
			for (size_t j = k; j < ncol; j++)
			{
				R[i][j] -= 2 * x[i] * x_R[j];
			}
		}

		// update Q
		for (size_t i = k; i < nrow; i++)
		{
			for (size_t j = 0; j < nrow; j++)
			{
				Q[i][j] += -2 * x[i] * x_Q[j];
			}
		}
	}


	// Transpose Q
	for (size_t i = 0; i < nrow; i++)
	{
		for (size_t j = i + 1; j < nrow; j++)
		{
			double tmp = Q[i][j];
			Q[i][j] = Q[j][i];
			Q[j][i] = tmp;
		}
	}

	/*
	// Check answer:
	for (size_t i = 0; i < nrow; i++)
	{
	for (size_t j = 0; j < ncol; j++)
	{
	double a_ij(0.0);
	for (size_t kk = 0; kk < nrow; kk++)
	{
	a_ij += Q[i][kk] * R[kk][j];
	}
	double err = fabs(a_ij - A[i][j]);
	if (err > 1E-4 * fabs(A[i][j]))
	{
	//vcm_cout << " High error in house holder QR method" << std::endl;
	}
	}
	}
	*/

	delete[] x; delete[] x_R; delete[] x_Q;
	return 0;
	#pragma endregion
}

// A = L LT is a square sparse matrix
bool MeshingLinearSolver::Cholesky_factorization(size_t num_rows, size_t** AJ, double** AV)
{	
	#pragma region Cholesky Factorization:
	for (size_t irow = 0; irow < num_rows; irow++)
	{
		double dot(0.0);
		if (irow > 0) dot = get_dot_product(AJ[irow], AV[irow], AJ[irow], AV[irow], 0, irow - 1);
		double aii(0.0);
		get_sparse_vector_entry(irow, aii, AJ[irow], AV[irow]);
		if (aii - dot < 1E-10)
		{
			return false;
		}
		double lii = sqrt(aii - dot);
		set_sparse_vector_entry(irow, lii, AJ[irow], AV[irow]);

		for (size_t jrow = irow + 1; jrow < num_rows; jrow++)
		{
			double dot(0.0);
			if (irow > 0) dot = get_dot_product(AJ[irow], AV[irow], AJ[jrow], AV[jrow], 0, irow - 1);
			double aji(0.0);
			get_sparse_vector_entry(irow, aji, AJ[jrow], AV[jrow]);
			double lji = (aji - dot) / lii;
			set_sparse_vector_entry(irow, lji, AJ[jrow], AV[jrow]);
		}
	}

	// removing the entries above the diagonal
	for (size_t irow = 0; irow < num_rows; irow++)
	{
		size_t nnz = AJ[irow][1];
		for (size_t i = 0; i < nnz; i++)
		{
			size_t col_index = AJ[irow][2 + i];
			if (col_index <= irow) continue;
			while (col_index > irow && i < nnz)
			{
				nnz--;
				AJ[irow][2 + i] = AJ[irow][2 + nnz];
				AV[irow][i] = AV[irow][nnz];
				col_index = AJ[irow][2 + i];
			}
		}
		AJ[irow][1] = nnz;
	}

	// Include L Transpose
	for (size_t irow = 1; irow < num_rows; irow++)
	{
		size_t nnz(AJ[irow][1]);
		for (size_t i = 0; i < nnz; i++)
		{
			size_t col_index = AJ[irow][2 + i];
			double val = AV[irow][i];
			if (col_index < irow)
			{
				set_sparse_vector_entry(irow, val, AJ[col_index], AV[col_index]);
			}
		}
	}
	return true;
	#pragma endregion
}

int MeshingLinearSolver::Cholesky_factorization_rank_one_update(size_t num_rows, size_t** CJ, double** CV, size_t* &xj, double* &xv)
{
	#pragma region Cholesky Factorization Rank one update:
	if (xj == 0) return 0;

	size_t nnz_x(xj[1]);
	size_t i(0);
	while (i < nnz_x)
	{
		bool advance_i(true);
		size_t irow = xj[2 + i];
		double xi_val = xv[i];

		double diag(0.0);
		get_sparse_vector_entry(irow, diag, CJ[irow], CV[irow]);
		double r = sqrt(diag * diag + xi_val * xi_val);
		double c = r / diag;
		double s = xi_val / diag;
		set_sparse_vector_entry(irow, r, CJ[irow], CV[irow]);

		size_t nnz(CJ[irow][1]);
		size_t pos = find_entry_position_binary(irow, CJ[irow], 2, 1 + nnz);
		while (CJ[irow][pos] <= irow && pos < 2 + nnz) pos++;
		size_t next_c_row(num_rows);
		if (pos < 2 + nnz) next_c_row = CJ[irow][pos];

		size_t xpos(3 + i);
		size_t next_x_row(num_rows);
		if (xpos < 2 + nnz_x) next_x_row = xj[xpos];

		while (next_c_row < num_rows || next_x_row < num_rows)
		{
			bool advance_pos(true), advance_xpos(true);
			double xj_val(0.0), lji_val(0.0);
			if (next_x_row <= next_c_row)
			{
				xj_val = xv[xpos - 2];
				if (next_c_row == next_x_row) lji_val = CV[irow][pos - 2];
				lji_val = (lji_val + s * xj_val) / c;
				set_sparse_vector_entry(irow, lji_val, CJ[next_x_row], CV[next_x_row]);
				set_sparse_vector_entry(next_x_row, lji_val, CJ[irow], CV[irow]);
				if (fabs(lji_val) < SPARSE_ENTRY_TOL) advance_pos = false;
				xj_val = c * xj_val - s * lji_val;
				xv[xpos - 2] = xj_val;
			}
			else
			{
				lji_val = CV[irow][pos - 2];
				lji_val = lji_val / c;
				set_sparse_vector_entry(irow, lji_val, CJ[next_c_row], CV[next_c_row]);
				set_sparse_vector_entry(next_c_row, lji_val, CJ[irow], CV[irow]);
				if (fabs(lji_val) < SPARSE_ENTRY_TOL) advance_pos = false;

				xj_val = -s * lji_val;
				if (fabs(xj_val) < SPARSE_ENTRY_TOL)
				{
					advance_xpos = false;
					if (next_c_row <= i) 
						advance_i = false;
				}
				set_sparse_vector_entry(next_c_row, xj_val, xj, xv);
			}
			
			if (advance_pos) pos++;

			nnz = CJ[irow][1];
			if (pos < 2 + nnz) next_c_row = CJ[irow][pos];
			else               next_c_row = num_rows;
			
			if (advance_xpos) xpos++;
			nnz_x = xj[1];
			if (xpos < 2 + nnz_x) next_x_row = xj[xpos];
			else                  next_x_row = num_rows;
		}
		if (advance_i) i++;
	}
	return 0;
	#pragma endregion
}

int MeshingLinearSolver::Cholesky_factorization_remove_row(size_t &num_rows, size_t** &CJ, double** &CV, size_t deleted_row_index)
{
	#pragma region Remove Row from a Cholseky factorization:

	if (deleted_row_index == num_rows - 1)
	{
		num_rows--;
		delete[] CJ[num_rows]; delete[] CV[num_rows];

		size_t** new_CJ = new size_t*[num_rows];
		double** new_CV = new double*[num_rows];
		for (size_t irow = 0; irow < num_rows; irow++)
		{
			new_CJ[irow] = 0; new_CV[irow] = 0;
			size_t nnz = CJ[irow][1];
			for (size_t j = 0; j < nnz; j++)
			{
				size_t col_index(CJ[irow][2 + j]);
				double val = CV[irow][j];
				if (col_index < deleted_row_index)
					set_sparse_vector_entry(col_index, val, new_CJ[irow], new_CV[irow]);
			}
			delete[] CJ[irow]; delete[] CV[irow];
		}
		delete[] CJ; delete[] CV;
		CJ = new_CJ; CV = new_CV;
		return 0;
	}

	size_t** S11_J(0); double** S11_V(0);
	if (deleted_row_index > 0)
	{
		S11_J = new size_t*[deleted_row_index];
		S11_V = new double*[deleted_row_index];
	}

	size_t** S13_J(0); double** S13_V(0);
	if (deleted_row_index > 0)
	{
		S13_J = new size_t*[deleted_row_index];
		S13_V = new double*[deleted_row_index];
	}

	// S11 and S13
	for (size_t irow = 0; irow < deleted_row_index; irow++)
	{
		S11_J[irow] = 0; S11_V[irow] = 0;
		S13_J[irow] = 0; S13_V[irow] = 0;
		size_t nnz(CJ[irow][1]);
		for (size_t j = 0; j < nnz; j++)
		{
			size_t col_index(CJ[irow][2 + j]);
			double val = CV[irow][j];
			if (col_index < deleted_row_index)
				set_sparse_vector_entry(col_index, val, S11_J[irow], S11_V[irow]);
			else if (col_index > deleted_row_index)
				set_sparse_vector_entry(col_index, val, S13_J[irow], S13_V[irow]);
		}
		delete[] CJ[irow]; delete[] CV[irow];
		CJ[irow] = 0; CV[irow] = 0;
	}

	size_t* s23_j(0); double* s23_v(0);
	if (true)
	{
		size_t nnz(CJ[deleted_row_index][1]);
		for (size_t j = 0; j < nnz; j++)
		{
			size_t col_index(CJ[deleted_row_index][2 + j]);
			double val = CV[deleted_row_index][j];
			if (col_index > deleted_row_index) set_sparse_vector_entry(col_index - deleted_row_index - 1, val, s23_j, s23_v);
		}
		delete[] CJ[deleted_row_index]; delete[] CV[deleted_row_index];
		CJ[deleted_row_index] = 0; CV[deleted_row_index] = 0;
	}

	size_t** S33_J(0); double** S33_V(0);
	if (deleted_row_index != num_rows - 1)
	{
		S33_J = new size_t*[num_rows - deleted_row_index - 1];
		S33_V = new double*[num_rows - deleted_row_index - 1];
	}

	for (size_t irow = deleted_row_index + 1; irow < num_rows; irow++)
	{
		S33_J[irow - deleted_row_index - 1] = 0; S33_V[irow - deleted_row_index - 1] = 0;
	}

	for (size_t irow = deleted_row_index + 1; irow < num_rows; irow++)
	{
		size_t nnz(CJ[irow][1]);
		for (size_t j = 0; j < nnz; j++)
		{
			size_t col_index(CJ[irow][2 + j]);
			double val = CV[irow][j];
			if (col_index > deleted_row_index)
			{
				set_sparse_vector_entry(col_index - deleted_row_index - 1, val, S33_J[irow - deleted_row_index - 1], S33_V[irow - deleted_row_index - 1]);
				set_sparse_vector_entry(irow - deleted_row_index - 1, val, S33_J[col_index - deleted_row_index - 1], S33_V[col_index - deleted_row_index - 1]);
			}
				
		}
		delete[] CJ[irow]; delete[] CV[irow];
		CJ[irow] = 0; CV[irow] = 0;
	}

	Cholesky_factorization_rank_one_update(num_rows - deleted_row_index - 1, S33_J, S33_V, s23_j, s23_v);

	delete[] s23_j; delete[] s23_v;

	for (size_t irow = 0; irow < deleted_row_index; irow++)
	{
		size_t nnz(S11_J[irow][1]);
		for (size_t j = 0; j < nnz; j++)
		{
			size_t col_index(S11_J[irow][2 + j]);
			double val = S11_V[irow][j];
			set_sparse_vector_entry(col_index, val, CJ[irow], CV[irow]);
		}
		delete[] S11_J[irow]; delete[] S11_V[irow];

		if (S13_J[irow] != 0)
		{
			nnz = S13_J[irow][1];
			for (size_t j = 0; j < nnz; j++)
			{
				size_t col_index(S13_J[irow][2 + j]);
				double val = S13_V[irow][j];
				set_sparse_vector_entry(col_index - 1, val, CJ[irow], CV[irow]);
				set_sparse_vector_entry(irow, val, CJ[col_index - 1], CV[col_index - 1]);
			}
			delete[] S13_J[irow]; delete[] S13_V[irow];
		}
	}
	if (S11_J != 0)
	{
		delete[] S11_J; delete[] S11_V;
	}
	
	if (S13_J != 0)
	{
		delete[] S13_J; delete[] S13_V;
	}

	for (size_t irow = deleted_row_index + 1; irow < num_rows; irow++)
	{
		size_t nnz(S33_J[irow - deleted_row_index - 1][1]);
		for (size_t j = 0; j < nnz; j++)
		{
			size_t col_index(S33_J[irow - deleted_row_index - 1][2 + j]);
			double val = S33_V[irow - deleted_row_index - 1][j];
			set_sparse_vector_entry(col_index + deleted_row_index, val, CJ[irow - 1], CV[irow - 1]);
		}
		delete[] S33_J[irow - deleted_row_index - 1]; delete[] S33_V[irow - deleted_row_index - 1];
	}

	if (S33_J != 0)
	{
		delete[] S33_J; delete[] S33_V;
	}
	
	num_rows--;
	return 0;
	#pragma endregion
}


bool MeshingLinearSolver::Cholesky_factorization_add_row(size_t& num_rows, size_t** &CJ, double** &CV, size_t* aj, double* av)
{
	#pragma region Add row to a Cholesky Facotrization: 
	// Solve L11 s_23 = a_12
	size_t* s12_j(0); double* s12_v(0);
	for (size_t i = 0; i < num_rows; i++)
	{
		double sum(0.0);
		if (s12_j != 0)
		{
			size_t nnz(s12_j[1]);
			for (size_t j = 0; j < nnz; j++)
			{
				size_t col_index(s12_j[2 + j]);
				if (col_index < i)
				{
					double val(0.0);
					get_sparse_vector_entry(col_index, val, CJ[i], CV[i]);
					sum += val * s12_v[j];
				}
			}
		}
		
		double bval(0.0);
		get_sparse_vector_entry(i, bval, aj, av);

		double diag(0.0);
		get_sparse_vector_entry(i, diag, CJ[i], CV[i]);

		set_sparse_vector_entry(i, (bval - sum) / diag, s12_j, s12_v);
	}

	double dot(0.0);
	if (s12_j != 0)
	{
		size_t nnz(s12_j[1]);
		for (size_t j = 0; j < nnz; j++)
		{
			dot += s12_v[j] * s12_v[j];
		}
	}
	
	double s_22(0.0);
	get_sparse_vector_entry(num_rows, s_22, aj, av);

	if (s_22 - dot < 1E-10)
	{
		delete[] s12_j; delete[] s12_v;
		//vcm_cout << "Warning:: Adding redundant row to Cholesky Factorization! diag = " << s_22 - dot << std::endl;
		return false;
	}

	s_22 = sqrt(s_22 - dot);

	size_t** CJ_new = new size_t*[num_rows + 1];
	double** CV_new = new double*[num_rows + 1];

	for (size_t i = 0; i < num_rows; i++)
	{
		CJ_new[i] = CJ[i];
		CV_new[i] = CV[i];
	}

	if (CJ != 0)
	{
		delete[] CJ; delete[] CV;
	}

	CJ_new[num_rows] = 0; CV_new[num_rows] = 0;

	// adding s12
	if (s12_j != 0)
	{
		size_t nnz(s12_j[1]);
		for (size_t j = 0; j < nnz; j++)
		{
			size_t row_index = s12_j[2 + j];
			double val = s12_v[j];
			set_sparse_vector_entry(num_rows, val, CJ_new[row_index], CV_new[row_index]);
			set_sparse_vector_entry(row_index, val, CJ_new[num_rows], CV_new[num_rows]);
		}
	}

	set_sparse_vector_entry(num_rows, s_22, CJ_new[num_rows], CV_new[num_rows]);

	if (s12_j != 0)
	{
		delete[] s12_j; delete[] s12_v;
	}
	
	CJ = CJ_new; CV = CV_new;
	num_rows++;
	return true;
	#pragma endregion
}


int MeshingLinearSolver::Cholesky_Solver(size_t num_rows, size_t** CJ, double** CV, double* b, double* x)
{
	#pragma region Cholesky Solver:
	double* y = new double[num_rows];
	// Solving L y = b
	for (size_t irow = 0; irow < num_rows; irow++)
	{
		double dot(0.0);
		if (irow > 0) dot = get_dot_product(CJ[irow], CV[irow], y, 0, irow - 1);
		double diag(0.0);
		get_sparse_vector_entry(irow, diag, CJ[irow], CV[irow]);
		y[irow] = (b[irow] - dot) / diag;
	}

	// Solving LT x = y
	size_t row_index(num_rows - 1);
	while (true)
	{
		double dot(0.0);
		dot = get_dot_product(CJ[row_index], CV[row_index], x, row_index + 1, num_rows - 1);
		double diag(0.0);
		get_sparse_vector_entry(row_index, diag, CJ[row_index], CV[row_index]);
		x[row_index] = (y[row_index] - dot) / diag;
		if (row_index == 0) break;
		row_index--;
	}
	delete[] y;
	return 0;
	#pragma endregion
}


int MeshingLinearSolver::set_sparse_vector_entry(size_t row, double val, size_t*& I, double*& V)
{
	#pragma region Set Sparse Matrix Entry:
	if (I == 0 && fabs(val) < SPARSE_ENTRY_TOL) return 1;

	if (I == 0)
	{
		I = new size_t[12];
		I[0] = 10; I[1] = 1; I[2] = row;
		V = new double[10]; V[0] = val;
		return 0;
	}

	if (I[1] == 0)
	{
		if (fabs(val) < SPARSE_ENTRY_TOL) return 0;

		I[1] = 1; I[2] = row; V[0] = val;
		return 0;
	}

	size_t cap(I[0]), num(I[1]);
	size_t pos = find_entry_position_binary(row, I, 2, 1 + num);

	if (I[pos] == row)
	{
		if (fabs(val) < SPARSE_ENTRY_TOL)
		{
			size_t last = 2 + num - 1;
			for (size_t j = pos; j < last; j++)
			{
				I[j] = I[j + 1];
				V[j - 2] = V[j - 1];
			}
			I[1]--;
		}
		else
		{
			V[pos - 2] = val;
		}
		return 0;
	}
	if (fabs(val) < SPARSE_ENTRY_TOL) return 1;

	if (num == cap)
	{
		size_t new_cap = cap * 2;
		size_t* vi = new size_t[new_cap + 2];
		double* v = new double[new_cap];
		vi[0] = new_cap; vi[1] = num;
		for (size_t i = 0; i < num; i++)
		{
			vi[i + 2] = I[i + 2];
			v[i] = V[i];
		}
		delete[] I; delete[] V;
		I = vi; V = v;
	}

	while (I[pos] > row && pos > 2) pos--;
	if (I[pos] < row) pos++;

	size_t last = 2 + num;
	for (size_t j = last; j > pos; j--)
	{
		I[j] = I[j - 1];
		V[j - 2] = V[j - 3];
	}
	I[pos] = row;
	V[pos - 2] = val;
	I[1]++;
	return 0;
	#pragma endregion
}

int MeshingLinearSolver::get_sparse_vector_entry(size_t row, double& val, size_t* I, double* V)
{
	#pragma region Get Sparse Vector Entry:
	val = 0.0;
	size_t num(I[1]);
	if (num == 0) return 0;

	size_t pos = find_entry_position_binary(row, I, 2, 1 + num);
	if (I[pos] == row)
		val = V[pos - 2];
	return 0; //  zero value
	#pragma endregion
}

size_t MeshingLinearSolver::find_entry_position_binary(size_t entry, size_t* I, size_t left, size_t right)
{
	#pragma region find using Binary search:

	size_t ipivot = (left + right) / 2;
	size_t pivot = I[ipivot];
	if (pivot == entry) return ipivot;
	if (left == right)  return ipivot;

	if (pivot < entry && ipivot + 1 <= right) return find_entry_position_binary(entry, I, ipivot + 1, right);
	if (pivot > entry && left + 1 <= ipivot) return find_entry_position_binary(entry, I, left, ipivot - 1);

	return ipivot;
	#pragma endregion
}

double MeshingLinearSolver::get_dot_product(size_t* I1, double* V1, size_t* I2, double* V2, size_t min_index, size_t max_index)
{
	#pragma region Dot product of Sparse vectors:
	double dot(0.0);

	if (I1 == 0 || I2 == 0) return 0.0;

	size_t nnz_1(I1[1]), nnz_2(I2[1]);
	if (nnz_1 < nnz_2)
	{
		for (size_t i = 0; i < nnz_1; i++)
		{
			size_t row_index(I1[2 + i]);
			if (row_index < min_index) continue;
			if (row_index > max_index) break;
			double val(0.0);
			get_sparse_vector_entry(row_index, val, I2, V2);
			dot += V1[i] * val;
		}
	}
	else
	{
		for (size_t i = 0; i < nnz_2; i++)
		{
			size_t row_index(I2[2 + i]);
			if (row_index < min_index) continue;
			if (row_index > max_index) break;
			double val(0.0);
			get_sparse_vector_entry(row_index, val, I1, V1);
			dot += val * V2[i];
		}
	}
	return dot;
	#pragma endregion
}

double MeshingLinearSolver::get_dot_product(size_t* I1, double* V1, double* V2, size_t min_index, size_t max_index)
{
	#pragma region Dot product of Sparse vectors:
	double dot(0.0);

	if (I1 == 0) return 0.0;

	size_t nnz_1(I1[1]);
	for (size_t i = 0; i < nnz_1; i++)
	{
		size_t row_index(I1[2 + i]);
		if (row_index < min_index) continue;
		if (row_index > max_index) break;
		dot += V1[i] * V2[row_index];
	}
	return dot;
	#pragma endregion
}
