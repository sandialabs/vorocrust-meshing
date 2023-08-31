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
//  MeshingRandomSampler.cpp                                      Last modified (02/22/2021) //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MeshingRandomSampler.h"

size_t MeshingRandomSampler::numOfRandomSamplers = 0;

MeshingRandomSampler::MeshingRandomSampler()
{
	numOfRandomSamplers++;
	size_t seed(time(0));
#ifdef FIXED_RNG_SEED
	seed = 1234567890;
#endif // FIXED_RNG_SEED
	if (numOfRandomSamplers > 0.001 * seed)
		seed += numOfRandomSamplers;
	else
		seed /= numOfRandomSamplers;

	// seed = 1614650530;

	//vcm_cout << "Random Seed = " << seed << std::endl;

	initiate_random_number_generator((unsigned long)(seed));
}


MeshingRandomSampler::MeshingRandomSampler(size_t thread_id)
{
	size_t seed(time(0));
#ifdef FIXED_RNG_SEED
	seed = 1234567890;
#endif // FIXED_RNG_SEED
	seed /= (thread_id + 1);
	initiate_random_number_generator((unsigned long)(seed));
}

MeshingRandomSampler::~MeshingRandomSampler()
{

}

void MeshingRandomSampler::initiate_random_number_generator(unsigned long x)
{
	#pragma region Initiate Random number generator:
    //assert(sizeof (double) >= 54) ;

    cc = 1.0 / 9007199254740992.0; // inverse of 2^53rd power
    size_t i;
    size_t qlen = indx = sizeof Q / sizeof Q[0];
    for (i = 0; i < qlen; i++) Q[i] = 0;

    c = 0.0;     /* current CSWB */
    zc = 0.0;	/* current SWB `borrow` */
    zx = 5212886298506819.0 / 9007199254740992.0;	/* SWB seed1 */
    zy = 2020898595989513.0 / 9007199254740992.0;	/* SWB seed2 */

    size_t j;
    double s, t;	 /* Choose 32 bits for x, 32 for y */
    if (x == 0) x = 123456789; /* default seeds */
    unsigned long y = 362436069; /* default seeds */

    /* Next, seed each Q[i], one bit at a time, */
    for (i = 0; i < qlen; i++)
    { /* using 9th bit from Cong+Xorshift */
        s = 0.0;
        t = 1.0;
        for (j = 0; j < 52; j++)
        {
            t = 0.5 * t; /* make t=.5/2^j */
            x = 69069 * x + 123;
            y ^= (y << 13);
            y ^= (y >> 17);
            y ^= (y << 5);
            if (((x + y) >> 23) & 1) s = s + t; /* change bit of s, maybe */
        }	 /* end j loop */
        Q[i] = s;
    } /* end i seed loop, Now generate 10^9 dUNI's: */
	#pragma endregion
}

double MeshingRandomSampler::generate_uniform_random_number()
{
	#pragma region Generate a Random Number:
    /* Takes 14 nanosecs, Intel Q6600,2.40GHz */
    int i, j;
    double t; /* t: first temp, then next CSWB value */
    /* First get zy as next lag-2 SWB */
    t = zx - zy - zc;
    zx = zy;
    if (t < 0)
    {
        zy = t + 1.0;
        zc = cc;
    }
    else
    {
        zy = t;
        zc = 0.0;
    }

    /* Then get t as the next lag-1220 CSWB value */
    if (indx < 1220)
        t = Q[indx++];
    else
    { /* refill Q[n] via Q[n-1220]-Q[n-1190]-c, */
        for (i = 0; i < 1220; i++)
        {
            j = (i < 30) ? i + 1190 : i - 30;
            t = Q[j] - Q[i] + c; /* Get next CSWB element */
            if (t > 0)
            {
                t = t - cc;
                c = cc;
            }
            else
            {
                t = t - cc + 1.0;
                c = 0.0;
            }
            Q[i] = t;
        }	 /* end i loop */
        indx = 1;
        t = Q[0]; /* set indx, exit 'else' with t=Q[0] */
    } /* end else segment; return t-zy mod 1 */

    return ((t < zy) ? 1.0 + (t - zy) : t - zy);
	#pragma endregion
}

int MeshingRandomSampler::sample_uniformly_from_box(size_t num_dim, double* xmin, double* xmax, double* dart)
{
	#pragma region Sample Uniformly From Box:
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		double u = generate_uniform_random_number();
		dart[idim] = xmin[idim] + u * (xmax[idim] - xmin[idim]);
	}
	return 0;
	#pragma endregion
}


int MeshingRandomSampler::sample_uniformly_from_unit_sphere(double* dart, size_t num_dim)
{
	#pragma region Sample a Point Unirofmly from a Sphere:
	size_t idim = 0;
	while (true)
	{
		double u1 = generate_uniform_random_number();
		double u2 = generate_uniform_random_number();
		double r = sqrt(-2 * log(u1));
		double theta = 2 * PI * u2;
		double n1 = r * cos(theta);
		double n2 = r * sin(theta);
		dart[idim] = n1; idim++; if (idim == num_dim) break;
		dart[idim] = n2; idim++; if (idim == num_dim) break;
	}
	double norm(0.0);
	for (idim = 0; idim < num_dim; idim++) norm += dart[idim] * dart[idim];

	norm = 1.0 / sqrt(norm);
	for (idim = 0; idim < num_dim; idim++) dart[idim] *= norm;
	return 0;
	#pragma endregion
}

int MeshingRandomSampler::sample_uniformly_from_simplex(double* dart, size_t num_dim, size_t num_corners, double** x)
{
	#pragma region Sample uniformly random from a simplex:
	double* u = new double[num_corners + 1];
	u[0] = 0.0;
	for (size_t idim = 1; idim < num_corners; idim++)
	{
		u[idim] = generate_uniform_random_number();
	}
	u[num_corners] = 1.0;

	quicksort(u, 0, num_corners);

	double* v = new double[num_corners];
	for (size_t ic = 0; ic < num_corners; ic++) v[ic] = (u[ic + 1] - u[ic]);

	for (size_t idim = 0; idim < num_dim; idim++) dart[idim] = 0.0;

	for (size_t ic = 0; ic < num_corners; ic++)
	{
		for (size_t idim = 0; idim < num_dim; idim++) dart[idim] += v[ic] * x[ic][idim];
	}

	delete[] u; delete[] v;

	return 0;
	#pragma endregion
}

size_t MeshingRandomSampler::sample_uniformly_from_discrete_cdf(size_t num_cells, double* cdf)
{
	#pragma region Sample From Discrete CDF:

	if (cdf[num_cells - 1] == 0.0)
	{
		//vcm_cout << "*** Warning zero cdf, num_cells = " << num_cells << ", cdf = " << cdf[num_cells - 1] << std::endl;
		size_t icell = size_t(num_cells * generate_uniform_random_number());
		if (icell == num_cells) icell--;
		return icell;
	}

	double u = generate_uniform_random_number() * cdf[num_cells - 1];
	size_t ist(0), iend(num_cells - 1);

	while (true)
	{
		#pragma region binary Serach to pick a cell:
		size_t imid = (ist + iend) / 2;
		double ulo(0.0);
		if (imid > 0) ulo = cdf[imid - 1];
		double uhi = cdf[imid];

		if (ulo <= u && u <= uhi)
		{
			return imid;
		}
		if (u < ulo) iend = imid;
		else         ist = imid + 1;
		#pragma endregion
	}
	return num_cells;
	#pragma endregion
}

int MeshingRandomSampler::quicksort(double* x, size_t left, size_t right)
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
