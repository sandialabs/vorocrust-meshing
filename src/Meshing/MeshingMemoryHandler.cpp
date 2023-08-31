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
// MeshingMemoryHandler.cpp                                       Last modified (07/12/2022) //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MeshingMemoryHandler.h"

MeshingMemoryHandler::MeshingMemoryHandler()
{	

}

MeshingMemoryHandler::~MeshingMemoryHandler()
{
   
}

int MeshingMemoryHandler::add_entry(size_t entry, size_t &num_entries, size_t* &I, size_t &capacity)
{
	#pragma region add an entry to a list:
	// add to neighbors
	I[num_entries] = entry;
	num_entries++;
	if (num_entries == capacity)
	{
		capacity *= 2;
		size_t* tmp = new size_t[capacity];
		for (size_t i = 0; i < num_entries; i++) tmp[i] = I[i];
		delete[] I;
		I = tmp;
	}
	return 0;
	#pragma endregion
}


int MeshingMemoryHandler::add_entry(bool* entry, size_t &num_entries, bool** &I, size_t &capacity)
{
	#pragma region add an entry to a list:
	// add to neighbors
	I[num_entries] = entry;
	num_entries++;
	if (num_entries == capacity)
	{
		capacity *= 2;
		bool** tmp = new bool*[capacity];
		for (size_t i = 0; i < num_entries; i++) tmp[i] = I[i];
		delete[] I;
		I = tmp;
	}
	return 0;
	#pragma endregion
}

int MeshingMemoryHandler::add_entry(size_t* entry, size_t &num_entries, size_t** &I, size_t &capacity)
{
	#pragma region add an entry to a list:
	// add to neighbors
	I[num_entries] = entry;
	num_entries++;
	if (num_entries == capacity)
	{
		capacity *= 2;
		size_t** tmp = new size_t*[capacity];
		for (size_t i = 0; i < num_entries; i++) tmp[i] = I[i];
		delete[] I;
		I = tmp;
	}
	return 0;
	#pragma endregion
}


int MeshingMemoryHandler::add_entry(double entry, size_t &num_entries, double* &I, size_t &capacity)
{
	#pragma region add an entry to a list:
	// add to neighbors
	I[num_entries] = entry;
	num_entries++;
	if (num_entries == capacity)
	{
		capacity *= 2;
		double* tmp = new double[capacity];
		for (size_t i = 0; i < num_entries; i++) tmp[i] = I[i];
		delete[] I;
		I = tmp;
	}
	return 0;
	#pragma endregion
}


int MeshingMemoryHandler::add_entry(double* entry, size_t &num_entries, double** &I, size_t &capacity)
{
	#pragma region add an entry to a list:
	// add to neighbors
	I[num_entries] = entry;
	num_entries++;
	if (num_entries == capacity)
	{
		capacity *= 2;
		double** tmp = new double*[capacity];
		for (size_t i = 0; i < num_entries; i++) tmp[i] = I[i];
		delete[] I;
		I = tmp;
	}
	return 0;
	#pragma endregion
}

size_t MeshingMemoryHandler::string_to_size_t(const std::string& s)
{
	std::istringstream i(s); size_t x;
	if (!(i >> x)) return 0; return x;
}

double MeshingMemoryHandler::string_to_double(const std::string& s)
{
	std::istringstream i(s); double x;
	if (!(i >> x)) return 0.0; return x;
}

int MeshingMemoryHandler::quicksort(size_t* x, size_t left, size_t right)
{
	#pragma region Quick Sort:
	size_t i = left, j = right;
	size_t pivot = x[(left + right) / 2];

	/* partition */
	while (i <= j)
	{
		while (x[i] < pivot)
			i++;
		while (x[j] > pivot)
			j--;
		if (i <= j)
		{
			size_t tmp = x[i];
			x[i] = x[j];
			x[j] = tmp;

			i++;
			if (j > 0) j--;
		}
	};

	/* recursion */

	if (j > 0 && left < j)
		quicksort(x, left, j);
	if (i < right)
		quicksort(x, i, right);

	return 0;
	#pragma endregion
}

int MeshingMemoryHandler::quicksort(double* x, size_t* I, size_t left, size_t right)
{
	#pragma region Quick Sort:
	size_t i = left, j = right;
	double pivot = x[(left + right) / 2];

	/* partition */
	while (i <= j)
	{
		while (x[i] < pivot)
			i++;
		while (x[j] > pivot)
			j--;
		if (i <= j)
		{
			double tmp = x[i];
			x[i] = x[j];
			x[j] = tmp;

			size_t tmpi = I[i];
			I[i] = I[j];
			I[j] = tmpi;

			i++;
			if (j > 0) j--;
		}
	};

	/* recursion */

	if (left < j && j > 0)
		quicksort(x, I, left, j);
	if (i < right)
		quicksort(x, I, i, right);

	return 0;
	#pragma endregion
}


int MeshingMemoryHandler::quicksort(size_t* x, size_t* I, size_t left, size_t right)
{
	#pragma region Quick Sort:
	size_t i = left, j = right;
	size_t pivot = x[(left + right) / 2];

	/* partition */
	while (i <= j)
	{
		while (x[i] < pivot)
			i++;
		while (x[j] > pivot)
			j--;
		if (i <= j)
		{
			size_t tmp = x[i];
			x[i] = x[j];
			x[j] = tmp;

			size_t tmpi = I[i];
			I[i] = I[j];
			I[j] = tmpi;

			i++;
			if (j > 0) j--;
		}
	};

	/* recursion */

	if (left < j && j > 0)
		quicksort(x, I, left, j);
	if (i < right)
		quicksort(x, I, i, right);

	return 0;
	#pragma endregion
}
