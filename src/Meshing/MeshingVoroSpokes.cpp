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
//  MeshingVoroSpokes.cpp                                         Last modified (05/15/2019) //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MeshingVoroSpokes.h"

MeshingVoroSpokes::MeshingVoroSpokes()
{
	_num_dim = 0; _num_seeds = 0; _budget = 0; _num_functions = 0; _desired_order = 0; _desired_orders = 0;
	_xmin = 0;  _xmax = 0; _sampling_function_index = 0;
	_seeds = 0; _f = 0; _cell_order = 0; _cell_weight = 0; _cell_cdf = 0; _cell_samples = 0; _cell_coef = 0;
	_vps_basis = Chebyshev; _use_total_order = false;
}

MeshingVoroSpokes::~MeshingVoroSpokes()
{
	clear_memory();
}

int MeshingVoroSpokes::init_vorospokes(size_t num_threads, size_t num_dim, size_t num_functions, size_t desired_order, double* xmin, double* xmax, bool adaptive_sampling_mode)
{
	#pragma region Init MeshingVoroSpokes:
	init_vorospokes(num_threads, num_dim, 100, num_functions, desired_order, xmin, xmax, adaptive_sampling_mode);
	return 0;
	#pragma endregion
}

int MeshingVoroSpokes::init_vorospokes(size_t num_threads, size_t num_dim, size_t budget, size_t num_functions, size_t desired_order, double* xmin, double* xmax, bool adaptive_sampling_mode)
{
	#pragma region Init MeshingVoroSpokes:
	_num_dim = num_dim; _budget = budget; _num_functions = num_functions; _desired_order = desired_order;
	_sampling_function_index = _num_functions;

	_xmin = new double[_num_dim];
	for (size_t idim = 0; idim < _num_dim; idim++) _xmin[idim] = xmin[idim];
	_xmax = new double[_num_dim];
	for (size_t idim = 0; idim < _num_dim; idim++) _xmax[idim] = xmax[idim];

	_seeds = new MeshingSmartTree(_num_dim);

	_f = new double*[_budget];
	for (size_t i = 0; i < _budget; i++) _f[i] = 0;
	_cell_coef = new double**[_budget];
	for (size_t i = 0; i < _budget; i++) _cell_coef[i] = 0;
	_cell_order = new size_t*[_budget];
	for (size_t i = 0; i < _budget; i++) _cell_order[i] = 0;

	_desired_orders = 0;


	if (adaptive_sampling_mode)
	{
		_cell_weight = new double[_budget];
		_cell_cdf = new double[_budget];
		_cell_samples = new double*[_budget];
		for (size_t i = 0; i < _budget; i++) _cell_samples[i] = 0;
	}
	else
	{
		_cell_weight = 0; _cell_cdf = 0; _cell_samples = 0;
	}

	init_random_samplers(num_threads);
	return 0;
	#pragma endregion
}

int MeshingVoroSpokes::clear_memory()
{
	#pragma region Clear Memory:
	if (_desired_orders != 0)
	{
		delete[] _desired_orders;
		_desired_orders = 0;
	}
	if (_xmin != 0)
	{
		delete[] _xmin;
		_xmin = 0;
	}
	if (_xmax != 0)
	{
		delete[] _xmax;
		_xmax = 0;
	}
	if (_seeds != 0)
	{
		delete _seeds;
		_seeds = 0;
	}

	if (_cell_weight != 0)
	{
		delete[] _cell_weight;
		_cell_weight = 0;
	}
	if (_cell_cdf != 0)
	{
		delete[] _cell_cdf;
		_cell_cdf = 0;
	}

	if (_f != 0)
	{
		for (size_t i = 0; i < _budget; i++)
		{
			if (_f[i] != 0) delete[] _f[i];
		}
		delete[] _f;
		_f = 0;
	}

	if (_cell_order != 0)
	{
		for (size_t i = 0; i < _budget; i++)
		{
			if (_cell_order[i] != 0) delete[] _cell_order[i];
		}
		delete[] _cell_order;
		_cell_order = 0;
	}

	if (_cell_samples != 0)
	{
		for (size_t i = 0; i < _budget; i++)
		{
			if (_cell_samples[i] != 0) delete[] _cell_samples[i];
		}
		delete[] _cell_samples;
		_cell_samples = 0;
	}

	if (_cell_coef != 0)
	{
		for (size_t i = 0; i < _budget; i++)
		{
			if (_cell_coef[i] != 0)
			{
				for (size_t j = 0; j < _num_functions; j++)
				{
					if (_cell_coef[i][j] != 0) delete[] _cell_coef[i][j];
				}
				delete[] _cell_coef[i];
			}
		}
		delete[] _cell_coef;
		_cell_coef = 0;
	}

	size_t num_threads = _rsamplers.size();
	for (size_t thread_id = 0; thread_id < num_threads; thread_id++)
	{
		delete _rsamplers[thread_id];
	}
	return 0;
	#pragma endregion
}

int MeshingVoroSpokes::set_sampling_function(size_t function_index, size_t num_spokes)
{
	#pragma region Set Sampling Function:
	_sampling_function_index = function_index;
	for (size_t iseed = 0; iseed < _num_seeds; iseed++)
	{
		integrate_vps(iseed, num_spokes, _sampling_function_index, _cell_weight[iseed], _cell_samples[iseed]);

		if (iseed == 0)
			_cell_cdf[iseed] = _cell_weight[_num_seeds];
		else
			_cell_cdf[iseed] = _cell_weight[_num_seeds] + _cell_cdf[iseed - 1];
	}
	return 0;
	#pragma endregion
}

int MeshingVoroSpokes::suggest_new_sample(double* xsample, size_t num_spokes)
{
	#pragma region Suggest New Seed:
	if (_num_seeds == 0)
	{
		_rsampler.sample_uniformly_from_box(_num_dim, _xmin, _xmax, xsample);
	}
	else
	{
		size_t selected_cell = _rsampler.sample_uniformly_from_discrete_cdf(_num_seeds, _cell_cdf);
		for (size_t idim = 0; idim < _num_dim; idim++) xsample[idim] = _cell_samples[selected_cell][idim];
		if (num_spokes > 0)
			integrate_vps(selected_cell, num_spokes, _sampling_function_index, _cell_weight[selected_cell], _cell_samples[selected_cell]);
	}
	return 0;
	#pragma endregion
}


int MeshingVoroSpokes::add_samples(size_t num_samples, double** samples_x, double** samples_f, size_t num_spokes)
{
	#pragma region Add Samples:
	for (size_t isample = 0; isample < num_samples; isample++)
	{
		_f[isample] = new double[_num_functions];
		for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) _f[isample][ifunc] = samples_f[isample][ifunc];
		_seeds->add_tree_point(_num_dim, samples_x[isample], 0, 0);
	}

	for (size_t isample = 0; isample < num_samples; isample++)
	{
		build_vps_model(isample, num_spokes);
	}
	_num_seeds = num_samples;
	return 0;
	#pragma endregion
}

// used in VoroCrust Sampler
int MeshingVoroSpokes::add_new_seed(double* xsample, double* fsample, size_t num_spokes)
{
	#pragma region Add New Seed:

	if (_cell_samples != 0 && _sampling_function_index != _num_functions)
	{
		// switch back to the adaptive sampling mode
		set_sampling_function(_num_functions, num_spokes);
	}

	_f[_num_seeds] = new double[_num_functions];
	for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) _f[_num_seeds][ifunc] = fsample[ifunc];

	if (_cell_samples != 0) _cell_samples[_num_seeds] = new double[_num_dim];

	if (_num_seeds == 0 || _desired_order == 0)
	{
		_seeds->add_tree_point(_num_dim, xsample, 0, 0);

		build_vps_model(_num_seeds, num_spokes);

		if (_cell_samples != 0) integrate_vps(_num_seeds, num_spokes, _sampling_function_index, _cell_weight[_num_seeds], _cell_samples[_num_seeds]);
		if (_cell_samples != 0) _cell_cdf[_num_seeds] = _cell_weight[_num_seeds];

		_num_seeds++;
	}
	else
	{
		size_t num_neighbors(0); size_t* neighbors(0); size_t* neighbor_hits(0);
		get_vps_neighbors(xsample, num_spokes, num_neighbors, neighbors, neighbor_hits);

		_seeds->add_tree_point(_num_dim, xsample, 0, 0);

		// update vps model of neighbor cells
		size_t num_threads = _rsamplers.size();
		size_t num_thread_neighbors = num_neighbors / num_threads;
		num_thread_neighbors++;
   	    bool last_seed_updated(false);

#if defined USE_OPEN_MP
		#pragma omp parallel
		{
			size_t thread_id = omp_get_thread_num();
#else
		for (size_t thread_id = 0; thread_id < num_threads; thread_id++)
		{
#endif
			for (size_t ineighbor = 0; ineighbor < num_thread_neighbors; ineighbor++)
			{

				size_t neighbor_index = thread_id * num_thread_neighbors + ineighbor;

				if (neighbor_index > num_neighbors) continue;

				size_t neighbor = _num_seeds;
				if (neighbor_index < num_neighbors)
				{
					neighbor = neighbors[neighbor_index];
				}
				if (neighbor_index == num_neighbors) last_seed_updated = true;

				build_vps_model(neighbor, num_spokes);
				if (_cell_samples != 0) integrate_vps(neighbor, num_spokes, _sampling_function_index, _cell_weight[neighbor], _cell_samples[neighbor]);
			}
		}

		if (!last_seed_updated)
		{
			build_vps_model(_num_seeds, num_spokes);
			if (_cell_samples != 0) integrate_vps(_num_seeds, num_spokes, _sampling_function_index, _cell_weight[_num_seeds], _cell_samples[_num_seeds]);
		}
		_num_seeds++;

		if (_cell_samples != 0)
		{
			// update cdf
			_cell_cdf[0] = _cell_weight[0];
			for (size_t i = 1; i < _num_seeds; i++) _cell_cdf[i] = _cell_cdf[i - 1] + _cell_weight[i];
		}
		delete[] neighbors; delete[] neighbor_hits;
	}
	if (_num_seeds == _budget) expand_vorospokes_containers();
	return 0;
	#pragma endregion
}

int MeshingVoroSpokes::evaluate_vps(double* x, double* fx)
{
	#pragma region Evaluate VPS:
	if (_num_seeds == 0)
	{
		for (size_t i = 0; i < _num_functions; i++) fx[i] = 0.0;
		return 0;
	}
	size_t iclosest; double hclosest(DBL_MAX);
	_seeds->get_closest_tree_point(x, 0, 0, iclosest, hclosest);
	return evaluate_vps(iclosest, x, fx);
	#pragma endregion
}

int MeshingVoroSpokes::evaluate_vps(double* x, double* fx, double &vps_err_est)
{
	#pragma region Evaluate VPS:
	if (_num_seeds == 0)
	{
		for (size_t i = 0; i < _num_functions; i++) fx[i] = 0.0;
		vps_err_est = DBL_MAX;
		return 0;
	}
	size_t iclosest; double hclosest(DBL_MAX);
	_seeds->get_closest_tree_point(x, 0, 0, iclosest, hclosest);
	if (_desired_orders != 0) return evaluate_vps(iclosest, x, fx);
	return evaluate_vps(iclosest, x, fx, vps_err_est);
	#pragma endregion
}

int MeshingVoroSpokes::evaluate_vwn(double* x, size_t num_spokes, size_t num_layers, double* fx, double &vwn_err_est)
{
	#pragma region Evaluate VWN:

	if (num_layers == 0) return evaluate_vps(x, fx, vwn_err_est); // plain VPS


	size_t num_neighbors(0); size_t* neighbors(0); size_t* neighbor_hits(0);

	_seeds->get_Voronoi_neighbors(num_spokes, num_layers, _num_dim, x, 0, 0, num_neighbors, neighbors, neighbor_hits);

	if (neighbor_hits[0] == SIZE_MAX)
	{
		size_t neighbor = neighbors[0];
		for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) fx[ifunc] = _f[neighbor][ifunc];
		delete[] neighbors; delete[] neighbor_hits;
		return 0;
	}

	for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) fx[ifunc] = 0.0;
	vwn_err_est = 0.0;

	double total_weight(0.0);
	double* ff = new double[_num_functions];
	for (size_t i = 0; i < num_neighbors; i++)
	{
		double angle = neighbor_hits[i] * 1.0 / num_spokes;
		size_t neighbor = neighbors[i];
		double* y = _seeds->get_tree_point(neighbor);

		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = x[idim] - y[idim];
			dst_sq += dx * dx;
		}

		double weight = angle * angle / dst_sq;

		double vps_err_estimate(0.0);

		evaluate_vps(neighbor, x, ff, vps_err_estimate);

		for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) fx[ifunc] += weight * ff[ifunc];
		vwn_err_est += weight * vps_err_estimate;

		total_weight += weight;
	}

	for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) fx[ifunc] /= total_weight;
	vwn_err_est /= total_weight;

	delete[] neighbors; delete[] neighbor_hits; delete[] ff;

	return 0;
	#pragma endregion
}

int MeshingVoroSpokes::plot_function(std::string file_name, size_t function_index, size_t num_contours)
{
	#pragma region Plot Solid Isocontours:
	//vcm_cout << ".: VPS Debug Mode :. Plotting ps files .... " << std::endl;

	std::fstream file(file_name.c_str(), std::ios::out);
	file << "%!PS-Adobe-3.0" << std::endl;
	file << "72 72 scale     % one unit = one inch" << std::endl;

	double xmin(_xmin[0]);
	double ymin(_xmin[1]);
	double Lx(_xmax[0] - _xmin[0]);
	double Ly(_xmax[1] - _xmin[1]);

	double scale_x, scale_y, scale;
	double shift_x, shift_y;

	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;

	if (scale_x < scale_y)
	{
		scale = scale_x;
		shift_x = 1.0 - xmin * scale;
		shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
	}
	else
	{
		scale = scale_y;
		shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
		shift_y = 1.0 - ymin * scale;
	}
	file << shift_x << " " << shift_y << " translate" << std::endl;

	file << "/Courier findfont" << std::endl;
	file << "0.12 scalefont" << std::endl;
	file << "setfont" << std::endl;

	std::vector<double> poly_x;
	std::vector<double> poly_y;

	size_t num_points = _seeds->get_num_tree_points();

	double f_min(DBL_MAX), f_max(-DBL_MAX);
	if (true)
	{
		double* dart = new double[2];
		double* dart_f = new double[_num_functions];
		size_t num_mc = 100000;
		for (size_t i = 0; i < num_mc; i++)
		{
			_rsampler.sample_uniformly_from_box(_num_dim, _xmin, _xmax, dart);
			double err_est(0.0);
			evaluate_vps(dart, dart_f, err_est);

			double fvps = err_est;
			if (function_index < _num_functions) fvps = dart_f[function_index];

			f_min = fmin(f_min, fvps);
			f_max = fmax(f_max, fvps);
		}
		delete[] dart; delete[] dart_f;
	}

	std::vector<double> contours;
	contours.push_back(f_min - 2 * (f_max - f_min));
	for (size_t i = 0; i < num_contours; i++) contours.push_back(f_min + (1.0 / num_contours) * i * (f_max - f_min));
	contours.push_back(f_max + 2 * (f_max - f_min));

	size_t n(500);
	double sx = Lx / n;
	double sy = Ly / n;

	for (size_t i = 0; i < n; i++)
	{
		double xo = _xmin[0] + i * sx;
		for (size_t j = 0; j < n; j++)
		{
			double fo(0.0), f1(0.0), f2(0.0), f3(0.0);

			double yo = _xmin[1] + j * sy;

			double* x = new double[2];
			double* fx = new double[_num_functions];

			x[0] = xo; x[1] = yo;
			evaluate_vps(x, fx, fo);
			if (function_index < _num_functions) fo = fx[function_index];

			x[0] += sx;
			evaluate_vps(x, fx, f1);
			if (function_index < _num_functions) f1 = fx[function_index];

			x[1] += sy;
			evaluate_vps(x, fx, f2);
			if (function_index < _num_functions) f2 = fx[function_index];

			x[0] -= sx;
			evaluate_vps(x, fx, f3);
			if (function_index < _num_functions) f3 = fx[function_index];;

			delete[] x; delete[] fx;

			size_t num_isocontours = contours.size();
			for (size_t icont = 1; icont < num_isocontours; icont++)
			{
				#pragma region Isocontouring:
				double contour = contours[icont];
				double contour_m = contours[icont - 1];

				//vcm_cout<< "contour_m = " << contour_m << " , contour = " << contour << std::endl;

				poly_x.clear(); poly_y.clear();

				// moving right
				if (fo >= contour_m - 1E-10 && fo < contour + 1E-10)
				{
					poly_x.push_back(xo);
					poly_y.push_back(yo);
					if ((fo > contour && f1 < contour) || (fo < contour && f1 > contour))
					{
						double h = sx * (contour - fo) / (f1 - fo);
						poly_x.push_back(xo + h);
						poly_y.push_back(yo);
					}
					else if ((fo > contour_m && f1 < contour_m) || (fo < contour_m && f1 > contour_m))
					{
						double h = sx * (contour_m - fo) / (f1 - fo);
						poly_x.push_back(xo + h);
						poly_y.push_back(yo);
					}
				}
				else if ((fo > contour_m && f1 < contour_m) || (fo < contour_m && f1 > contour_m))
				{
					double hm = sx * (contour_m - fo) / (f1 - fo);
					double h = hm;
					if ((fo > contour && f1 < contour) || (fo < contour && f1 > contour))
					{
						h = sx * (contour - fo) / (f1 - fo);
					}
					if (h < hm)
					{
						double tmp = h; h = hm; hm = tmp;
					}
					poly_x.push_back(xo + hm);
					poly_y.push_back(yo);

					if (h - hm > 1E-10)
					{
						poly_x.push_back(xo + h);
						poly_y.push_back(yo);
					}
				}
				else if ((fo > contour && f1 < contour) || (fo < contour && f1 > contour))
				{
					double h = sx * (contour - fo) / (f1 - fo);
					poly_x.push_back(xo + h);
					poly_y.push_back(yo);
				}

				// moving up
				if (f1 >= contour_m - 1E-10 && f1 < contour + 1E-10)
				{
					poly_x.push_back(xo + sx);
					poly_y.push_back(yo);
					if ((f1 > contour && f2 < contour) || (f1 < contour && f2 > contour))
					{
						double h = sy * (contour - f1) / (f2 - f1);
						poly_x.push_back(xo + sx);
						poly_y.push_back(yo + h);
					}
					else if ((f1 > contour_m && f2 < contour_m) || (f1 < contour_m && f2 > contour_m))
					{
						double h = sy * (contour_m - f1) / (f2 - f1);
						poly_x.push_back(xo + sx);
						poly_y.push_back(yo + h);
					}

				}
				else if ((f1 > contour_m && f2 < contour_m) || (f1 < contour_m && f2 > contour_m))
				{
					double hm = sy * (contour_m - f1) / (f2 - f1);
					double h = hm;
					if ((f1 > contour && f2 < contour) || (f1 < contour && f2 > contour))
					{
						h = sy * (contour - f1) / (f2 - f1);
					}
					if (h < hm)
					{
						double tmp = h; h = hm; hm = tmp;
					}
					poly_x.push_back(xo + sx);
					poly_y.push_back(yo + hm);

					if (h - hm > 1E-10)
					{
						poly_x.push_back(xo + sx);
						poly_y.push_back(yo + h);
					}
				}
				else if ((f1 > contour && f2 < contour) || (f1 < contour && f2 > contour))
				{
					double h = sy * (contour - f1) / (f2 - f1);
					poly_x.push_back(xo + sx);
					poly_y.push_back(yo + h);
				}

				// moving left
				if (f2 >= contour_m - 1E-10 && f2 < contour + 1E-10)
				{
					poly_x.push_back(xo + sx);
					poly_y.push_back(yo + sy);
					if ((f2 > contour && f3 < contour) || (f2 < contour && f3 > contour))
					{
						double h = sx * (contour - f2) / (f3 - f2);
						poly_x.push_back(xo + sx - h);
						poly_y.push_back(yo + sy);
					}
					else if ((f2 > contour_m && f3 < contour_m) || (f2 < contour_m && f3 > contour_m))
					{
						double h = sx * (contour_m - f2) / (f3 - f2);
						poly_x.push_back(xo + sx - h);
						poly_y.push_back(yo + sy);
					}
				}
				else if ((f2 > contour_m && f3 < contour_m) || (f2 < contour_m && f3 > contour_m))
				{
					double hm = sx * (contour_m - f2) / (f3 - f2);
					double h = hm;
					if ((f2 > contour && f3 < contour) || (f2 < contour && f3 > contour))
					{
						h = sx * (contour - f2) / (f3 - f2);
					}
					if (h < hm)
					{
						double tmp = h; h = hm; hm = tmp;
					}
					poly_x.push_back(xo + sx - hm);
					poly_y.push_back(yo + sy);

					if (h - hm > 1E-10)
					{
						poly_x.push_back(xo + sx - h);
						poly_y.push_back(yo + sy);
					}
				}
				else if ((f2 > contour && f3 < contour) || (f2 < contour && f3 > contour))
				{
					double h = sx * (contour - f2) / (f3 - f2);
					poly_x.push_back(xo + sx - h);
					poly_y.push_back(yo + sy);
				}

				// moving down
				if (f3 >= contour_m - 1E-10 && f3 < contour + 1E-10)
				{
					poly_x.push_back(xo);
					poly_y.push_back(yo + sy);
					if ((f3 > contour && fo < contour) || (f3 < contour && fo > contour))
					{
						double h = sy * (contour - f3) / (fo - f3);
						poly_x.push_back(xo);
						poly_y.push_back(yo + sy - h);
					}
					else if ((f3 > contour_m && fo < contour_m) || (f3 < contour_m && fo > contour_m))
					{
						double h = sy * (contour_m - f3) / (fo - f3);
						poly_x.push_back(xo);
						poly_y.push_back(yo + sy - h);
					}
				}
				else if ((f3 > contour_m && fo < contour_m) || (f3 < contour_m && fo > contour_m))
				{
					double hm = sy * (contour_m - f3) / (fo - f3);
					double h = hm;
					if ((f3 > contour && fo < contour) || (f3 < contour && fo > contour))
					{
						h = sy * (contour - f3) / (fo - f3);
					}
					if (h < hm)
					{
						double tmp = h; h = hm; hm = tmp;
					}
					poly_x.push_back(xo);
					poly_y.push_back(yo + sy - hm);

					if (h - hm > 1E-10)
					{
						poly_x.push_back(xo);
						poly_y.push_back(yo + sy - h);
					}
				}
				else if ((f3 > contour && fo < contour) || (f3 < contour && fo > contour))
				{
					double h = sy * (contour - f3) / (fo - f3);
					poly_x.push_back(xo);
					poly_y.push_back(yo + sy - h);
				}


				size_t num_corners(poly_x.size());
				if (num_corners > 1)
				{
					double gs = 1.0 - icont * 1.0 / num_isocontours;
					file << "newpath" << std::endl;
					file << poly_x[0] * scale << " " << poly_y[0] * scale << " moveto" << std::endl;
					//vcm_cout<< "*** x = " <<  poly_x[0] << ", y = " << poly_y[0] << std::endl;
					for (size_t icorner = 1; icorner < num_corners; icorner++)
					{
						file << poly_x[icorner] * scale << " " << poly_y[icorner] * scale << " lineto" << std::endl;
						//vcm_cout << "*** x = " <<  poly_x[icorner] << ", y = " << poly_y[icorner] << std::endl;
					}
					//vcm_cout << std::endl;

					file << "closepath" << std::endl;
					file << "gsave" << std::endl;
					file << "grestore" << std::endl;

					double r, g, b;

					if (gs < 0.25)     r = 1.0;
					//else if (gs < 0.5) r = 2.0 - 4.0 * gs;
					else if (gs < 0.5) r = 1.0 - 16.0 * (gs - 0.25) * (gs - 0.25);
					else               r = 0.0;

					double go(0.25), gn(1.0 - go);
					if (gs < go)      g = gs / go;
					else if (gs < gn) g = 1.0;
					else              g = 1.0 / (1.0 - gn) - gs / (1.0 - gn);


					if (gs < 0.5)       b = 0.0;
					else if (gs < 0.75) b = 1.0 - 16.0 * (gs - 0.75) * (gs - 0.75);
					else                b = 1.0;

					file << r << " " << g << " " << b << " setrgbcolor" << std::endl;

					file << " fill" << std::endl;
				}
				#pragma endregion
			}
		}
	}

	for (size_t i = 0; i < num_points; i++)
	{
		#pragma region Plot Points:
		double* x = _seeds->get_tree_point(i);

		file << "newpath" << std::endl;
		file << x[0] * scale << " " << x[1] * scale << " " << 0.0005 * scale << " 0 360 arc" << std::endl;
		file << "closepath" << std::endl;
		file << "gsave" << std::endl;
		file << "0 0 0 setrgbcolor" << std::endl; // Inactive seeds
		file << "fill" << std::endl;
		file << "grestore" << std::endl;
		file << "0.0 setlinewidth" << std::endl;
		file << "0 0 0 setrgbcolor" << std::endl; // discs borders
		file << "stroke " << std::endl;

		std::stringstream ss;
		ss << i;
		std::string str = ss.str();

		//file << "newpath " << x[i][0] * scale << " " << x[i][1] * scale << " moveto (" << str << ") show" << std::endl;
		#pragma endregion
	}
	file << "showpage" << std::endl;
	return 0;
	#pragma endregion
}



int MeshingVoroSpokes::plot_function_slice(std::string file_name, double* po, double slice_width, size_t dim_1, size_t dim_2, size_t function_index, size_t num_contours, double f_min, double f_max, bool use_vwn, bool plot_points)
{
	#pragma region Plot Solid Isocontours:
	//vcm_cout << ".: VPS Debug Mode :. Plotting ps files .... " << std::endl;

	std::fstream file(file_name.c_str(), std::ios::out);
	file << "%!PS-Adobe-3.0" << std::endl;
	file << "72 72 scale     % one unit = one inch" << std::endl;

	double xmin(_xmin[dim_1]);
	double ymin(_xmin[dim_2]);
	double Lx(_xmax[dim_1] - _xmin[dim_1]);
	double Ly(_xmax[dim_2] - _xmin[dim_2]);

	double scale_x, scale_y, scale;
	double shift_x, shift_y;

	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;

	if (scale_x < scale_y)
	{
		scale = scale_x;
		shift_x = 1.0 - xmin * scale;
		shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
	}
	else
	{
		scale = scale_y;
		shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
		shift_y = 1.0 - ymin * scale;
	}
	file << shift_x << " " << shift_y << " translate" << std::endl;

	file << "/Courier findfont" << std::endl;
	file << "0.12 scalefont" << std::endl;
	file << "setfont" << std::endl;

	std::vector<double> poly_x;
	std::vector<double> poly_y;

	size_t num_points = _seeds->get_num_tree_points();

	std::vector<double> contours;
	contours.push_back(f_min - 2 * (f_max - f_min));
	for (size_t i = 0; i < num_contours; i++) contours.push_back(f_min + (1.0 / num_contours) * i * (f_max - f_min));
	contours.push_back(f_max + 2 * (f_max - f_min));

	size_t n(100);
	double sx = Lx / n;
	double sy = Ly / n;

	double err_est(0.0);
	for (size_t i = 0; i < n; i++)
	{
		double xo = _xmin[dim_1] + i * sx;
		for (size_t j = 0; j < n; j++)
		{
			double fo(0.0), f1(0.0), f2(0.0), f3(0.0);

			double yo = _xmin[dim_2] + j * sy;

			double* x = new double[_num_dim];
			for (size_t idim = 0; idim < _num_dim; idim++) x[idim] = po[idim];

			double* fx = new double[_num_functions];

			x[dim_1] = xo; x[dim_2] = yo;
			if (use_vwn) evaluate_vwn(x, 100, 1, fx, err_est);
			else         evaluate_vps(x, fx, err_est);
			fo = fx[function_index];
			if (function_index < _num_functions) fo = fx[function_index];

			x[dim_1] += sx;
			if (use_vwn) evaluate_vwn(x, 100, 1, fx, err_est);
			else         evaluate_vps(x, fx, err_est);
			if (function_index < _num_functions) f1 = fx[function_index];

			x[dim_2] += sy;
			if (use_vwn) evaluate_vwn(x, 100, 1, fx, err_est);
			else         evaluate_vps(x, fx, err_est);
			if (function_index < _num_functions) f2 = fx[function_index];

			x[dim_1] -= sx;
			if (use_vwn) evaluate_vwn(x, 100, 1, fx, err_est);
			else         evaluate_vps(x, fx, err_est);
			if (function_index < _num_functions) f3 = fx[function_index];;

			delete[] x; delete[] fx;

			size_t num_isocontours = contours.size();
			for (size_t icont = 1; icont < num_isocontours; icont++)
			{
				#pragma region Isocontouring:
				double contour = contours[icont];
				double contour_m = contours[icont - 1];

				//vcm_cout<< "contour_m = " << contour_m << " , contour = " << contour << std::endl;

				poly_x.clear(); poly_y.clear();

				// moving right
				if (fo >= contour_m - 1E-10 && fo < contour + 1E-10)
				{
					poly_x.push_back(xo);
					poly_y.push_back(yo);
					if ((fo > contour && f1 < contour) || (fo < contour && f1 > contour))
					{
						double h = sx * (contour - fo) / (f1 - fo);
						poly_x.push_back(xo + h);
						poly_y.push_back(yo);
					}
					else if ((fo > contour_m && f1 < contour_m) || (fo < contour_m && f1 > contour_m))
					{
						double h = sx * (contour_m - fo) / (f1 - fo);
						poly_x.push_back(xo + h);
						poly_y.push_back(yo);
					}
				}
				else if ((fo > contour_m && f1 < contour_m) || (fo < contour_m && f1 > contour_m))
				{
					double hm = sx * (contour_m - fo) / (f1 - fo);
					double h = hm;
					if ((fo > contour && f1 < contour) || (fo < contour && f1 > contour))
					{
						h = sx * (contour - fo) / (f1 - fo);
					}
					if (h < hm)
					{
						double tmp = h; h = hm; hm = tmp;
					}
					poly_x.push_back(xo + hm);
					poly_y.push_back(yo);

					if (h - hm > 1E-10)
					{
						poly_x.push_back(xo + h);
						poly_y.push_back(yo);
					}
				}
				else if ((fo > contour && f1 < contour) || (fo < contour && f1 > contour))
				{
					double h = sx * (contour - fo) / (f1 - fo);
					poly_x.push_back(xo + h);
					poly_y.push_back(yo);
				}

				// moving up
				if (f1 >= contour_m - 1E-10 && f1 < contour + 1E-10)
				{
					poly_x.push_back(xo + sx);
					poly_y.push_back(yo);
					if ((f1 > contour && f2 < contour) || (f1 < contour && f2 > contour))
					{
						double h = sy * (contour - f1) / (f2 - f1);
						poly_x.push_back(xo + sx);
						poly_y.push_back(yo + h);
					}
					else if ((f1 > contour_m && f2 < contour_m) || (f1 < contour_m && f2 > contour_m))
					{
						double h = sy * (contour_m - f1) / (f2 - f1);
						poly_x.push_back(xo + sx);
						poly_y.push_back(yo + h);
					}

				}
				else if ((f1 > contour_m && f2 < contour_m) || (f1 < contour_m && f2 > contour_m))
				{
					double hm = sy * (contour_m - f1) / (f2 - f1);
					double h = hm;
					if ((f1 > contour && f2 < contour) || (f1 < contour && f2 > contour))
					{
						h = sy * (contour - f1) / (f2 - f1);
					}
					if (h < hm)
					{
						double tmp = h; h = hm; hm = tmp;
					}
					poly_x.push_back(xo + sx);
					poly_y.push_back(yo + hm);

					if (h - hm > 1E-10)
					{
						poly_x.push_back(xo + sx);
						poly_y.push_back(yo + h);
					}
				}
				else if ((f1 > contour && f2 < contour) || (f1 < contour && f2 > contour))
				{
					double h = sy * (contour - f1) / (f2 - f1);
					poly_x.push_back(xo + sx);
					poly_y.push_back(yo + h);
				}

				// moving left
				if (f2 >= contour_m - 1E-10 && f2 < contour + 1E-10)
				{
					poly_x.push_back(xo + sx);
					poly_y.push_back(yo + sy);
					if ((f2 > contour && f3 < contour) || (f2 < contour && f3 > contour))
					{
						double h = sx * (contour - f2) / (f3 - f2);
						poly_x.push_back(xo + sx - h);
						poly_y.push_back(yo + sy);
					}
					else if ((f2 > contour_m && f3 < contour_m) || (f2 < contour_m && f3 > contour_m))
					{
						double h = sx * (contour_m - f2) / (f3 - f2);
						poly_x.push_back(xo + sx - h);
						poly_y.push_back(yo + sy);
					}
				}
				else if ((f2 > contour_m && f3 < contour_m) || (f2 < contour_m && f3 > contour_m))
				{
					double hm = sx * (contour_m - f2) / (f3 - f2);
					double h = hm;
					if ((f2 > contour && f3 < contour) || (f2 < contour && f3 > contour))
					{
						h = sx * (contour - f2) / (f3 - f2);
					}
					if (h < hm)
					{
						double tmp = h; h = hm; hm = tmp;
					}
					poly_x.push_back(xo + sx - hm);
					poly_y.push_back(yo + sy);

					if (h - hm > 1E-10)
					{
						poly_x.push_back(xo + sx - h);
						poly_y.push_back(yo + sy);
					}
				}
				else if ((f2 > contour && f3 < contour) || (f2 < contour && f3 > contour))
				{
					double h = sx * (contour - f2) / (f3 - f2);
					poly_x.push_back(xo + sx - h);
					poly_y.push_back(yo + sy);
				}

				// moving down
				if (f3 >= contour_m - 1E-10 && f3 < contour + 1E-10)
				{
					poly_x.push_back(xo);
					poly_y.push_back(yo + sy);
					if ((f3 > contour && fo < contour) || (f3 < contour && fo > contour))
					{
						double h = sy * (contour - f3) / (fo - f3);
						poly_x.push_back(xo);
						poly_y.push_back(yo + sy - h);
					}
					else if ((f3 > contour_m && fo < contour_m) || (f3 < contour_m && fo > contour_m))
					{
						double h = sy * (contour_m - f3) / (fo - f3);
						poly_x.push_back(xo);
						poly_y.push_back(yo + sy - h);
					}
				}
				else if ((f3 > contour_m && fo < contour_m) || (f3 < contour_m && fo > contour_m))
				{
					double hm = sy * (contour_m - f3) / (fo - f3);
					double h = hm;
					if ((f3 > contour && fo < contour) || (f3 < contour && fo > contour))
					{
						h = sy * (contour - f3) / (fo - f3);
					}
					if (h < hm)
					{
						double tmp = h; h = hm; hm = tmp;
					}
					poly_x.push_back(xo);
					poly_y.push_back(yo + sy - hm);

					if (h - hm > 1E-10)
					{
						poly_x.push_back(xo);
						poly_y.push_back(yo + sy - h);
					}
				}
				else if ((f3 > contour && fo < contour) || (f3 < contour && fo > contour))
				{
					double h = sy * (contour - f3) / (fo - f3);
					poly_x.push_back(xo);
					poly_y.push_back(yo + sy - h);
				}


				size_t num_corners(poly_x.size());
				if (num_corners > 1)
				{
					double gs = 1.0 - icont * 1.0 / num_isocontours;
					file << "newpath" << std::endl;
					file << poly_x[0] * scale << " " << poly_y[0] * scale << " moveto" << std::endl;
					//vcm_cout<< "*** x = " <<  poly_x[0] << ", y = " << poly_y[0] << std::endl;
					for (size_t icorner = 1; icorner < num_corners; icorner++)
					{
						file << poly_x[icorner] * scale << " " << poly_y[icorner] * scale << " lineto" << std::endl;
						//vcm_cout << "*** x = " <<  poly_x[icorner] << ", y = " << poly_y[icorner] << std::endl;
					}
					//vcm_cout << std::endl;

					file << "closepath" << std::endl;
					file << "gsave" << std::endl;
					file << "grestore" << std::endl;

					double r, g, b;

					if (gs < 0.25)     r = 1.0;
					//else if (gs < 0.5) r = 2.0 - 4.0 * gs;
					else if (gs < 0.5) r = 1.0 - 16.0 * (gs - 0.25) * (gs - 0.25);
					else               r = 0.0;

					double go(0.25), gn(1.0 - go);
					if (gs < go)      g = gs / go;
					else if (gs < gn) g = 1.0;
					else              g = 1.0 / (1.0 - gn) - gs / (1.0 - gn);


					if (gs < 0.5)       b = 0.0;
					else if (gs < 0.75) b = 1.0 - 16.0 * (gs - 0.75) * (gs - 0.75);
					else                b = 1.0;

					file << r << " " << g << " " << b << " setrgbcolor" << std::endl;

					file << " fill" << std::endl;
				}
				#pragma endregion
			}
		}
	}
	if (plot_points)
	{
		#pragma region Plot Points:
		for (size_t i = 0; i < num_points; i++)
		{
			double* x = _seeds->get_tree_point(i);

			double hproj(0.0);
			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				if (idim == dim_1 || idim == dim_2) continue;
				double dx = x[idim] - po[idim];
				hproj += dx * dx;
			}

			if (hproj < 0.25 * slice_width * slice_width)
			{
				file << "newpath" << std::endl;
				file << x[dim_1] * scale << " " << x[dim_2] * scale << " " << 0.0005 * scale << " 0 360 arc" << std::endl;
				file << "closepath" << std::endl;
				file << "gsave" << std::endl;
				file << "0 0 0 setrgbcolor" << std::endl; // Inactive seeds
				file << "fill" << std::endl;
				file << "grestore" << std::endl;
				file << "0.0 setlinewidth" << std::endl;
				file << "0 0 0 setrgbcolor" << std::endl; // discs borders
				file << "stroke " << std::endl;

				//std::stringstream ss;
				//ss << i;
				//std::string str = ss.str();
				//file << "newpath " << x[dim_1] * scale << " " << x[dim_2] * scale << " moveto (" << str << ") show" << std::endl;
			}
		}
		#pragma endregion
	}

	file << "showpage" << std::endl;
	return 0;
	#pragma endregion
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Private Methods:
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int MeshingVoroSpokes::init_random_samplers(size_t num_threads)
{
	#pragma region Init Random Samplers:
	_rsamplers.resize(num_threads);
	for (size_t thread_id = 0; thread_id < num_threads; thread_id++)
	{
		_rsamplers[thread_id] = new MeshingRandomSampler(thread_id);
	}
	_seeds->init_random_samplers(num_threads);
	return 0;
	#pragma endregion
}

int MeshingVoroSpokes::expand_vorospokes_containers()
{
	#pragma region Expand Vorospokes Containers:
	size_t new_budget = _budget * 2;
	double** f = new double*[new_budget];
	double*** cell_coef = new double**[new_budget];
	size_t** cell_order = new size_t * [new_budget];

	double* cell_weight(0);
	double* cell_cdf(0);
	double** cell_samples(0);

	if (_cell_weight != 0) cell_weight = new double[new_budget];
	if (_cell_cdf != 0)    cell_cdf = new double[new_budget];
	if (_cell_samples != 0) cell_samples = new double* [new_budget];

	for (size_t i = 0; i < _budget; i++) f[i] = _f[i];
	for (size_t i = 0; i < _budget; i++) cell_coef[i] = _cell_coef[i];
	for (size_t i = 0; i < _budget; i++) cell_order[i] = _cell_order[i];

	for (size_t i = _budget; i < new_budget; i++) f[i] = 0;
	for (size_t i = _budget; i < new_budget; i++) cell_coef[i] = 0;
	for (size_t i = _budget; i < new_budget; i++) cell_order[i] = 0;

	delete[] _f; delete[] _cell_coef;  delete[] _cell_order;

	if (_cell_weight != 0) { for (size_t i = 0; i < _budget; i++) cell_weight[i] = _cell_weight[i]; }
	if (_cell_cdf != 0) { for (size_t i = 0; i < _budget; i++) cell_cdf[i] = _cell_cdf[i]; }
	if (_cell_samples != 0) { for (size_t i = 0; i < _budget; i++) cell_samples[i] = _cell_samples[i]; }

	if (_cell_weight != 0) delete[] _cell_weight;
	if (_cell_cdf != 0) delete[] _cell_cdf;
	if (_cell_samples != 0) delete[] _cell_samples;

	if (_cell_samples != 0) { for (size_t i = _budget; i < new_budget; i++) cell_samples[i] = 0; }

	_budget = new_budget;
	_cell_weight = cell_weight; _cell_cdf = cell_cdf; _f = f;
	_cell_order = cell_order; _cell_samples = cell_samples; _cell_coef = cell_coef;

	return 0;
	#pragma endregion
}

double MeshingVoroSpokes::evaluate_basis(double x, size_t ibasis)
{
	if (_vps_basis == monomials) return pow(x, ibasis);
	if (_vps_basis == Chebyshev)
	{
		if (ibasis == 0) return 1.0;
		else if (ibasis == 1) return x;

		double fo = 1.0; double f1 = x; size_t jbasis(2);
		return evaluate_basis_chebyshev(x, f1, fo, ibasis, jbasis);
	}
	return 0;
}

double MeshingVoroSpokes::evaluate_basis_chebyshev(double x, double fm, double fmm, size_t ibasis, size_t& jbasis)
{
	double f = 2 * x * fm - fmm;
	if (ibasis == jbasis) return f;
	jbasis++; fmm = fm; fm = f;
	return evaluate_basis_chebyshev(x, fm, fmm, ibasis, jbasis);
}

size_t MeshingVoroSpokes::get_num_basis(size_t num_dim, size_t order)
{
	#pragma region Get Number of Basis:

	size_t * index = new size_t[num_dim];
	for (size_t idim = 0; idim < num_dim; idim++) index[idim] = 0;

	size_t num_basis(0);
	while (true)
	{
		size_t kdim = num_dim - 1;
		index[kdim]++;

		while (kdim > 0 && index[kdim] > order)
		{
			index[kdim] = 0;
			index[kdim - 1]++;
			kdim--;
		}
		if (index[0] > order) break;

		size_t sum(0);
		for (size_t idim = 0; idim < num_dim; idim++) sum += index[idim];

		if (!_use_total_order && sum > order) continue;

		num_basis++;
	}
	delete[] index;
	return num_basis;
	#pragma endregion
}


double MeshingVoroSpokes::evaluate_basis(size_t ibasis, size_t num_dim, size_t order, double* x, size_t &basis_order)
{
	#pragma region Evaluate Basis:
	size_t * index = new size_t[num_dim];
	for (size_t idim = 0; idim < num_dim; idim++) index[idim] = 0;

	size_t num_basis(0);
	basis_order = 0;
	while (true)
	{
		size_t kdim = num_dim - 1;
		index[kdim]++;

		while (kdim > 0 && index[kdim] > order)
		{
			index[kdim] = 0;
			index[kdim - 1]++;
			kdim--;
		}
		if (index[0] > order) break;

		size_t sum(0);
		for (size_t idim = 0; idim < num_dim; idim++) sum += index[idim];

		if (!_use_total_order && sum > order) continue;

		if (num_basis == ibasis)
		{
			double prod(1.0);
			for (size_t idim = 0; idim < num_dim; idim++)
			{
				prod *= evaluate_basis(x[idim], index[idim]);
			}
			for (size_t idim = 0; idim < num_dim; idim++) basis_order += index[idim];
			delete[] index;
			return prod;
		}
		num_basis++;
	}
	delete[] index;
	return 0.0;
	#pragma endregion
}

int MeshingVoroSpokes::evaluate_vps(size_t icell, double* x, double* f_x)
{
	#pragma region Evaluate VPS:

	double* y = _seeds->get_tree_point(icell);

	for (size_t ifunc = 0; ifunc < _num_functions; ifunc++)
	{
		double f_vps = _f[icell][ifunc]; // vps zero order using neighbor cell

		double min_err_coef(1E-4);
		if (f_vps == DBL_MAX) min_err_coef = 1E-2; // boost the error of a singular cell

		size_t num_basis = get_num_basis(_num_dim, _cell_order[icell][ifunc]);
		double* xx = new double[_num_dim];
		for (size_t idim = 0; idim < _num_dim; idim++) xx[idim] = x[idim] - y[idim];
		for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
		{
			size_t basis_order(0);
			double fbasis = evaluate_basis(ibasis, _num_dim, _cell_order[icell][ifunc], xx, basis_order);
			f_vps += (_cell_coef[icell][ifunc][ibasis] * fbasis);
		}
		delete[] xx;

		f_x[ifunc] = f_vps;
	}
	return 0;
	#pragma endregion
}

int MeshingVoroSpokes::evaluate_vps(size_t icell, double* x, double* f_x, double &vps_err_est)
{
	#pragma region Evaluate VPS:
	vps_err_est = 0.0;
	double* y = _seeds->get_tree_point(icell);

	for (size_t ifunc = 0; ifunc < _num_functions; ifunc++)
	{
		double f_vps = _f[icell][ifunc]; // vps zero order using neighbor cell

		double min_err_coef(1E-4);
		if (f_vps == DBL_MAX) min_err_coef = 1E-2; // boost the error of a singular cell

		if (_cell_order[icell][ifunc] == 0)
		{
			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				double dx = y[idim] - x[idim];
				vps_err_est += min_err_coef * fabs(dx);
			}
		}
		else
		{
			size_t num_basis = get_num_basis(_num_dim, _cell_order[icell][ifunc]);
			double* xx = new double[_num_dim];
			for (size_t idim = 0; idim < _num_dim; idim++) xx[idim] = x[idim] - y[idim];
			for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
			{
				size_t basis_order;
				double fbasis = evaluate_basis(ibasis, _num_dim, _cell_order[icell][ifunc], xx, basis_order);
				if (basis_order == _cell_order[icell][ifunc])
				{
					double err_coef = fabs(_cell_coef[icell][ifunc][ibasis]);
					if (err_coef < min_err_coef) err_coef = min_err_coef;
					vps_err_est += fabs(err_coef * fbasis);
				}
				f_vps += (_cell_coef[icell][ifunc][ibasis] * fbasis);
			}
			delete[] xx;
		}
		f_x[ifunc] = f_vps;
	}
	return 0;
	#pragma endregion
}

int MeshingVoroSpokes::build_vps_model(size_t icell, size_t num_spokes)
{
	#pragma region Build VPS Model via weighted regression:

	if (_cell_order[icell] == 0) _cell_order[icell] = new size_t[_num_functions];

	if (_cell_coef[icell] == 0)
	{
		_cell_coef[icell] = new double*[_num_functions];
		for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) _cell_coef[icell][ifunc] = 0;
	}

	for (size_t ifunc = 0; ifunc < _num_functions; ifunc++)
	{
		if (_cell_coef[icell][ifunc] != 0) delete[] _cell_coef[icell][ifunc];
		_cell_coef[icell][ifunc] = 0;

		if (_desired_order == 0)
		{
			_cell_order[icell][ifunc] = 0;
			continue;
		}

		size_t num_neighbors(0); size_t* neighbors(0); size_t* neighbor_hits(0);
		double* x = _seeds->get_tree_point(icell);
		get_vps_neighbors(x, num_spokes, num_neighbors, neighbors, neighbor_hits);

		size_t num_basis = get_num_basis(_num_dim, _desired_order);
		_cell_order[icell][ifunc] = _desired_order;
		while (num_neighbors < 2 * num_basis && _cell_order[icell][ifunc] > 0)
		{
			_cell_order[icell][ifunc]--;
			num_basis = get_num_basis(_num_dim, _cell_order[icell][ifunc]);
		}

		if (_cell_order[icell][ifunc] == 0)
		{
			delete[] neighbors; delete[] neighbor_hits;
			continue;
		}

		_cell_coef[icell][ifunc] = new double[num_basis];
		for (size_t ibasis = 0; ibasis < num_basis; ibasis++) _cell_coef[icell][ifunc][ibasis] = 0.0;


		double** A = new double*[num_neighbors - 1]; double* b = new double[num_neighbors - 1];
		double* xx = new double[_num_dim];

		for (size_t i = 1; i < num_neighbors; i++)
		{
			A[i - 1] = new double[num_basis];

			size_t neighbor = neighbors[i];
			double* y = _seeds->get_tree_point(neighbor);

			for (size_t idim = 0; idim < _num_dim; idim++) xx[idim] = y[idim] - x[idim];

			double angle = neighbor_hits[i] * 1.0 / num_spokes;

			double dst_sq(0.0);
			for (size_t idim = 0; idim < _num_dim; idim++) dst_sq += xx[idim] * xx[idim];

			if (fabs(dst_sq < 1E-12))
			{
				vcm_cout << "MeshingVoroSpokes_WARNING:: Vonoronoi seeds " << icell << " and " << neighbor << " are too close to each other, distance = " << sqrt(dst_sq) << std::endl;
			}

			double weight = angle * angle / dst_sq;

			for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
			{
				size_t basis_order(0);
				A[i - 1][ibasis] = weight * evaluate_basis(ibasis, _num_dim, _cell_order[icell][ifunc], xx, basis_order);
			}
			b[i - 1] = weight * (_f[neighbor][ifunc] - _f[icell][ifunc]);
		}

		// CG_LS(num_neighbors - 1, num_basis, A, b, 1000, 1E-6, _cell_coef[icell][ifunc]);

		_lin.LS_QR_Solver(num_neighbors - 1, num_basis, A, b, _cell_coef[icell][ifunc]);

		for (size_t ipoint = 0; ipoint < num_neighbors - 1; ipoint++) delete[] A[ipoint];
		delete[] A; delete[] b; delete[] xx;

		delete[] neighbors; delete[] neighbor_hits;
	}
	return 0;
	#pragma endregion
}

int MeshingVoroSpokes::build_vps_model(size_t icell, size_t num_basis, double** c)
{
	#pragma region Build VPS Model via weighted regression:
	if (_cell_order[icell] == 0) _cell_order[icell] = new size_t[_num_functions];

	if (_cell_coef[icell] == 0)
	{
		_cell_coef[icell] = new double*[_num_functions];
		for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) _cell_coef[icell][ifunc] = 0;
	}

	for (size_t ifunc = 0; ifunc < _num_functions; ifunc++)
	{
		if (_cell_coef[icell][ifunc] != 0) delete[] _cell_coef[icell][ifunc];
		_cell_coef[icell][ifunc] = 0;
		if (_desired_order == 0 && _desired_orders == 0)
		{
			_cell_order[icell][ifunc] = 0;
			continue;
		}
		_cell_order[icell][ifunc] = _desired_order;
		_cell_coef[icell][ifunc] = new double[num_basis];
		for (size_t ibasis = 0; ibasis < num_basis; ibasis++) _cell_coef[icell][ifunc][ibasis] = c[ifunc][ibasis];
	}
	return 0;
	#pragma endregion
}

int MeshingVoroSpokes::get_vps_neighbors(double* x, size_t num_spokes, size_t &num_neighbors, size_t* &neighbors, size_t* &neighbor_hits)
{
	#pragma region Get VPS Neighbors:
	size_t num_layers(1);
	size_t num_basis = get_num_basis(_num_dim, _desired_order);
	while (true)
	{
		size_t num_old_neighbors(num_neighbors);
		_seeds->get_Voronoi_neighbors(num_spokes, num_layers, _num_dim, x, 0, 0, num_neighbors, neighbors, neighbor_hits);
		if (num_old_neighbors == num_neighbors) break; //no new neighbors
		if (num_neighbors >= 2 * num_basis) break;
		delete[] neighbors; delete[] neighbor_hits;
		neighbors = 0; neighbor_hits = 0;
		num_layers++;
	}
	return 0;
	#pragma endregion
}

double MeshingVoroSpokes::spherical_integration_constant(size_t num_dim)
{
	#pragma region Sepherical Integration Constant:
	if (num_dim / 2 * 2 == num_dim)
	{
		// num_dim is an even number
		size_t k = num_dim / 2;
		size_t gamma(k);
		while (k > 1)
		{
			k--; gamma *= k;
		}
		return num_dim * pow(PI, num_dim / 2) / gamma;
	}
	else
	{
		// num_dim is an odd number
		size_t k(num_dim); double gamma(0.5 * k);
		while (k > 2)
		{
			k -= 2; gamma *= (0.5 * k);
		}
		return num_dim * pow(PI, num_dim / 2) / gamma;
	}
	return 0.0;
	#pragma endregion
}

int MeshingVoroSpokes::integrate_vps(size_t icell, size_t num_spokes, size_t function_index, double &integration_val, double* xsample)
{
	#pragma region Integrate VPS and sample a point from the Underlying distribution:

	double* x = _seeds->get_tree_point(icell);

	size_t num_neighbors; size_t* neighbors; size_t* neighbor_hits;
	double** spoke_end_points = new double*[num_spokes];
	size_t* spoke_neighbors = new size_t[num_spokes];

	double* ff = new double[_num_functions];

	_seeds->get_Voronoi_neighbors(num_spokes, 1, _num_dim, x, spoke_end_points, spoke_neighbors, num_neighbors, neighbors, neighbor_hits);

	double* spoke_weights = new double[num_spokes];
	for (size_t ispoke = 0; ispoke < num_spokes; ispoke++)
	{
		#pragma region Estimate Spoke Weights:
		double spoke_length(0.0);

		if (spoke_neighbors[ispoke] == SIZE_MAX)
		{
			double umin = DBL_MAX;
			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				if (spoke_end_points[ispoke][idim] > 1E-10)
				{
					double u = (_xmax[idim] - x[idim]) / spoke_end_points[ispoke][idim];
					if (u < umin) umin = u;
				}
				else if (spoke_end_points[ispoke][idim] < -1E-10)
				{
					double u = (_xmin[idim] - x[idim]) / spoke_end_points[ispoke][idim];
					if (u < umin) umin = u;
				}
			}
			spoke_length = umin;

			// adjust spoke end_points
			for (size_t idim = 0; idim < _num_dim; idim++) spoke_end_points[ispoke][idim] = x[idim] + umin * spoke_end_points[ispoke][idim];
		}
		else
		{
			// ensure that spoke end is inside the bounding box
			double umin(1.0);
			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				if (spoke_end_points[ispoke][idim] < _xmin[idim])
				{
					double u = (x[idim] - _xmin[idim]) / (x[idim] - spoke_end_points[ispoke][idim]);
					if (u < umin) umin = u;
				}
				if (spoke_end_points[ispoke][idim] > _xmax[idim])
				{
					double u = (_xmax[idim] - x[idim]) / (spoke_end_points[ispoke][idim] - x[idim]);
					if (u < umin) umin = u;
				}
			}

			if (umin < 1.0)
			{
				for (size_t idim = 0; idim < _num_dim; idim++)
				{
					double dx = spoke_end_points[ispoke][idim] - x[idim];
					spoke_end_points[ispoke][idim] = x[idim] + umin * dx;
				}
			}

			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				double dx = spoke_end_points[ispoke][idim] - x[idim];
				spoke_length += dx * dx;
			}
			spoke_length = sqrt(spoke_length);
		}

		if (true)
		{
			size_t spoke_order(_desired_order);
			if (function_index < _num_functions) spoke_order = _cell_order[icell][function_index];
			else
			{
				spoke_order = 0;
				for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) spoke_order = std::max(spoke_order, _cell_order[icell][ifunc]);
			}

			if (spoke_order == 0) spoke_order++;

			size_t num_points(spoke_order + 1), num_spoke_basis(num_points);

			double* y = new double[_num_dim];
			double** A = new double*[num_points];
			double* b = new double[num_points];
			double dr = spoke_length / (num_points - 1);
			for (size_t ipoint = 0; ipoint < num_points; ipoint++)
			{
				double r = ipoint * dr;
				for (size_t idim = 0; idim < _num_dim; idim++) y[idim] = x[idim] + (r / spoke_length) * (spoke_end_points[ispoke][idim] - x[idim]);

				A[ipoint] = new double[spoke_order + 1];
				for (size_t ibasis = 0; ibasis < num_spoke_basis; ibasis++)
				{
					A[ipoint][ibasis] = pow(r, ibasis);
				}

				double vps_err_est(0.0);
				evaluate_vps(icell, y, ff, vps_err_est);

				if (function_index < _num_functions)
					b[ipoint] = ff[function_index];
				else
					b[ipoint] = vps_err_est;
			}

			double* cr = new double[num_spoke_basis];
			_lin.LS_QR_Solver(num_points, num_spoke_basis, A, b, cr);

			double val(0.0);
			for (size_t ibasis = 0; ibasis < num_spoke_basis; ibasis++)
			{
				size_t power = ibasis + _num_dim;
				val += cr[ibasis] * pow(spoke_length, power) / power;
				if (function_index == _num_functions) val -= cr[ibasis] * pow(0.5 * spoke_length, power) / power;
			}
			val *= spherical_integration_constant(_num_dim);

			spoke_weights[ispoke] = val;

			for (size_t ipoint = 0; ipoint < num_points; ipoint++) delete[] A[ipoint];
			delete[] A; delete[] b; delete[] y; delete[] cr;
		}
		#pragma endregion
	}

	integration_val = 0.0;
	for (size_t ispoke = 0; ispoke < num_spokes; ispoke++)
	{
		if (spoke_weights[ispoke] < 0.0) spoke_weights[ispoke] = 0.0;
		integration_val += spoke_weights[ispoke];
	}
	integration_val /= num_spokes;

	if (xsample != 0)
	{
		#pragma region Sample a point from the underlying distribution:
		for (size_t ispoke = 1; ispoke < num_spokes; ispoke++) spoke_weights[ispoke] += spoke_weights[ispoke - 1];

		size_t selected_spoke = _rsampler.sample_uniformly_from_discrete_cdf(num_spokes, spoke_weights);

		// pick a point along the chosen spoke
		double u = _rsampler.generate_uniform_random_number();
		if (integration_val < 1E-10)
		{
			// sample uniformly
			double t = pow(u, 1.0 / _num_dim);
			if (function_index == _num_functions && t < 0.5) t = 0.5;
			for (size_t idim = 0; idim < _num_dim; idim++) xsample[idim] = x[idim] + t * (spoke_end_points[selected_spoke][idim] - x[idim]);
		}
		else
		{
			double spoke_length(0.0);
			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				double dx = spoke_end_points[selected_spoke][idim] - x[idim];
				spoke_length += dx * dx;
			}
			spoke_length = sqrt(spoke_length);

			size_t spoke_order(_desired_order);
			if (function_index < _num_functions) spoke_order = _cell_order[icell][function_index];
			else
			{
				spoke_order = 0;
				for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) spoke_order = std::max(spoke_order, _cell_order[icell][ifunc]);
			}

			if (spoke_order == 0) spoke_order++;

			size_t num_points(spoke_order + 1), num_spoke_basis(num_points);

			double* y = new double[_num_dim];
			double** A = new double*[num_points];
			double* b = new double[num_points];
			double dr = spoke_length / (num_points - 1);
			for (size_t ipoint = 0; ipoint < num_points; ipoint++)
			{
				double r = ipoint * dr;
				for (size_t idim = 0; idim < _num_dim; idim++) y[idim] = x[idim] + (r / spoke_length) * (spoke_end_points[selected_spoke][idim] - x[idim]);

				A[ipoint] = new double[spoke_order + 1];
				for (size_t ibasis = 0; ibasis < num_spoke_basis; ibasis++)
				{
					A[ipoint][ibasis] = pow(r, ibasis);
				}

				double vps_err_est(0.0);
				evaluate_vps(icell, y, ff, vps_err_est);

				if (function_index < _num_functions)
					b[ipoint] = ff[function_index];
				else
					b[ipoint] = vps_err_est;
			}

			double* cr = new double[num_spoke_basis];
			_lin.LS_QR_Solver(num_points, num_spoke_basis, A, b, cr);

			u = _rsampler.generate_uniform_random_number();

			double val(0.0), val_max(0.0);
			for (size_t ibasis = 0; ibasis < num_spoke_basis; ibasis++)
			{
				size_t power = ibasis + _num_dim;
				val_max += cr[ibasis] * pow(spoke_length, power) / power;
				if (function_index == _num_functions) val_max -= cr[ibasis] * pow(0.5 * spoke_length, power) / power;
			}

			double tlo(0.0), thi(1.0);
			if (function_index == _num_functions) tlo = 0.5;
			while (true)
			{
				#pragma region Binary search for the point along the spoke that corresponds to u:
				double t = 0.5 * (tlo + thi);
				double val(0.0);
				for (size_t ibasis = 0; ibasis < num_spoke_basis; ibasis++)
				{
					size_t power = ibasis + _num_dim;
					val += cr[ibasis] * pow(t * spoke_length, power) / power;
					if (function_index == _num_functions) val -= cr[ibasis] * pow(0.5 * spoke_length, power) / power;
				}
				val /= val_max;

				if (fabs(val - u) < 1E-6)
				{
					for (size_t idim = 0; idim < _num_dim; idim++)
					{
						xsample[idim] = x[idim] + t * (spoke_end_points[selected_spoke][idim] - x[idim]);
					}
					break;
				}

				if (u < val) thi = t;
				else         tlo = t;
				#pragma endregion
			}
			for (size_t ipoint = 0; ipoint < num_points; ipoint++) delete[] A[ipoint];
			delete[] A; delete[] b; delete[] y; delete[] cr;
		}
		#pragma endregion
	}

	for (size_t ispoke = 0; ispoke < num_spokes; ispoke++) delete[] spoke_end_points[ispoke];
	delete[] spoke_end_points; delete[] spoke_neighbors;  delete[] spoke_weights;
	delete[] neighbors; delete[] neighbor_hits;
	delete[] ff;

	return 0;
	#pragma endregion
}


int MeshingVoroSpokes::integrate_vps(size_t icell, size_t num_spokes, size_t function_index, double& integration_val, double* xsample, size_t thread_id)
{
	#pragma region Integrate VPS and sample a point from the Underlying distribution:

	double* x = _seeds->get_tree_point(icell);

	size_t num_neighbors; size_t* neighbors; size_t* neighbor_hits;
	double** spoke_end_points = new double* [num_spokes];
	size_t* spoke_neighbors = new size_t[num_spokes];

	double* ff = new double[_num_functions];

	_seeds->get_Voronoi_neighbors(thread_id, num_spokes, 1, _num_dim, x, spoke_end_points, spoke_neighbors, num_neighbors, neighbors, neighbor_hits);

	double* spoke_weights = new double[num_spokes];
	for (size_t ispoke = 0; ispoke < num_spokes; ispoke++)
	{
		#pragma region Estimate Spoke Weights:
		double spoke_length(0.0);

		if (spoke_neighbors[ispoke] == SIZE_MAX)
		{
			double umin = DBL_MAX;
			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				if (spoke_end_points[ispoke][idim] > 1E-10)
				{
					double u = (_xmax[idim] - x[idim]) / spoke_end_points[ispoke][idim];
					if (u < umin) umin = u;
				}
				else if (spoke_end_points[ispoke][idim] < -1E-10)
				{
					double u = (_xmin[idim] - x[idim]) / spoke_end_points[ispoke][idim];
					if (u < umin) umin = u;
				}
			}
			spoke_length = umin;

			// adjust spoke end_points
			for (size_t idim = 0; idim < _num_dim; idim++) spoke_end_points[ispoke][idim] = x[idim] + umin * spoke_end_points[ispoke][idim];
		}
		else
		{
			// ensure that spoke end is inside the bounding box
			double umin(1.0);
			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				if (spoke_end_points[ispoke][idim] < _xmin[idim])
				{
					double u = (x[idim] - _xmin[idim]) / (x[idim] - spoke_end_points[ispoke][idim]);
					if (u < umin) umin = u;
				}
				if (spoke_end_points[ispoke][idim] > _xmax[idim])
				{
					double u = (_xmax[idim] - x[idim]) / (spoke_end_points[ispoke][idim] - x[idim]);
					if (u < umin) umin = u;
				}
			}

			if (umin < 1.0)
			{
				for (size_t idim = 0; idim < _num_dim; idim++)
				{
					double dx = spoke_end_points[ispoke][idim] - x[idim];
					spoke_end_points[ispoke][idim] = x[idim] + umin * dx;
				}
			}

			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				double dx = spoke_end_points[ispoke][idim] - x[idim];
				spoke_length += dx * dx;
			}
			spoke_length = sqrt(spoke_length);
		}

		if (true)
		{
			size_t spoke_order(_desired_order);
			if (function_index < _num_functions) spoke_order = _cell_order[icell][function_index];
			else
			{
				spoke_order = 0;
				for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) spoke_order = std::max(spoke_order, _cell_order[icell][ifunc]);
			}

			if (spoke_order == 0) spoke_order++;

			size_t num_points(spoke_order + 1), num_spoke_basis(num_points);

			double* y = new double[_num_dim];
			double** A = new double* [num_points];
			double* b = new double[num_points];
			double dr = spoke_length / (num_points - 1);
			for (size_t ipoint = 0; ipoint < num_points; ipoint++)
			{
				double r = ipoint * dr;
				for (size_t idim = 0; idim < _num_dim; idim++) y[idim] = x[idim] + (r / spoke_length) * (spoke_end_points[ispoke][idim] - x[idim]);

				A[ipoint] = new double[spoke_order + 1];
				for (size_t ibasis = 0; ibasis < num_spoke_basis; ibasis++)
				{
					A[ipoint][ibasis] = pow(r, ibasis);
				}

				double vps_err_est(0.0);
				evaluate_vps(icell, y, ff, vps_err_est);

				if (function_index < _num_functions)
					b[ipoint] = ff[function_index];
				else
					b[ipoint] = vps_err_est;
			}

			double* cr = new double[num_spoke_basis];
			_lin.LS_QR_Solver(num_points, num_spoke_basis, A, b, cr);

			double val(0.0);
			for (size_t ibasis = 0; ibasis < num_spoke_basis; ibasis++)
			{
				size_t power = ibasis + _num_dim;
				val += cr[ibasis] * pow(spoke_length, power) / power;
				if (function_index == _num_functions) val -= cr[ibasis] * pow(0.5 * spoke_length, power) / power;
			}
			val *= spherical_integration_constant(_num_dim);

			spoke_weights[ispoke] = val;

			for (size_t ipoint = 0; ipoint < num_points; ipoint++) delete[] A[ipoint];
			delete[] A; delete[] b; delete[] y; delete[] cr;
		}
		#pragma endregion
	}

	integration_val = 0.0;
	for (size_t ispoke = 0; ispoke < num_spokes; ispoke++)
	{
		if (spoke_weights[ispoke] < 0.0) spoke_weights[ispoke] = 0.0;
		integration_val += spoke_weights[ispoke];
	}
	integration_val /= num_spokes;

	if (xsample != 0)
	{
		#pragma region Sample a point from the underlying distribution:
		for (size_t ispoke = 1; ispoke < num_spokes; ispoke++) spoke_weights[ispoke] += spoke_weights[ispoke - 1];

		size_t selected_spoke = _rsamplers[thread_id]->sample_uniformly_from_discrete_cdf(num_spokes, spoke_weights);

		// pick a point along the chosen spoke
		double u = _rsamplers[thread_id]->generate_uniform_random_number();
		if (integration_val < 1E-10)
		{
			// sample uniformly
			double t = pow(u, 1.0 / _num_dim);
			if (function_index == _num_functions && t < 0.5) t = 0.5;
			for (size_t idim = 0; idim < _num_dim; idim++) xsample[idim] = x[idim] + t * (spoke_end_points[selected_spoke][idim] - x[idim]);
		}
		else
		{
			double spoke_length(0.0);
			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				double dx = spoke_end_points[selected_spoke][idim] - x[idim];
				spoke_length += dx * dx;
			}
			spoke_length = sqrt(spoke_length);

			size_t spoke_order(_desired_order);
			if (function_index < _num_functions) spoke_order = _cell_order[icell][function_index];
			else
			{
				spoke_order = 0;
				for (size_t ifunc = 0; ifunc < _num_functions; ifunc++) spoke_order = std::max(spoke_order, _cell_order[icell][ifunc]);
			}

			if (spoke_order == 0) spoke_order++;

			size_t num_points(spoke_order + 1), num_spoke_basis(num_points);

			double* y = new double[_num_dim];
			double** A = new double* [num_points];
			double* b = new double[num_points];
			double dr = spoke_length / (num_points - 1);
			for (size_t ipoint = 0; ipoint < num_points; ipoint++)
			{
				double r = ipoint * dr;
				for (size_t idim = 0; idim < _num_dim; idim++) y[idim] = x[idim] + (r / spoke_length) * (spoke_end_points[selected_spoke][idim] - x[idim]);

				A[ipoint] = new double[spoke_order + 1];
				for (size_t ibasis = 0; ibasis < num_spoke_basis; ibasis++)
				{
					A[ipoint][ibasis] = pow(r, ibasis);
				}

				double vps_err_est(0.0);
				evaluate_vps(icell, y, ff, vps_err_est);

				if (function_index < _num_functions)
					b[ipoint] = ff[function_index];
				else
					b[ipoint] = vps_err_est;
			}

			double* cr = new double[num_spoke_basis];
			_lin.LS_QR_Solver(num_points, num_spoke_basis, A, b, cr);

			u = _rsamplers[thread_id]->generate_uniform_random_number();

			double val(0.0), val_max(0.0);
			for (size_t ibasis = 0; ibasis < num_spoke_basis; ibasis++)
			{
				size_t power = ibasis + _num_dim;
				val_max += cr[ibasis] * pow(spoke_length, power) / power;
				if (function_index == _num_functions) val_max -= cr[ibasis] * pow(0.5 * spoke_length, power) / power;
			}

			double tlo(0.0), thi(1.0);
			if (function_index == _num_functions) tlo = 0.5;
			while (true)
			{
				#pragma region Binary search for the point along the spoke that corresponds to u:
				double t = 0.5 * (tlo + thi);
				double val(0.0);
				for (size_t ibasis = 0; ibasis < num_spoke_basis; ibasis++)
				{
					size_t power = ibasis + _num_dim;
					val += cr[ibasis] * pow(t * spoke_length, power) / power;
					if (function_index == _num_functions) val -= cr[ibasis] * pow(0.5 * spoke_length, power) / power;
				}
				val /= val_max;

				if (fabs(val - u) < 1E-6)
				{
					for (size_t idim = 0; idim < _num_dim; idim++)
					{
						xsample[idim] = x[idim] + t * (spoke_end_points[selected_spoke][idim] - x[idim]);
					}
					break;
				}

				if (u < val) thi = t;
				else         tlo = t;
				#pragma endregion
			}
			for (size_t ipoint = 0; ipoint < num_points; ipoint++) delete[] A[ipoint];
			delete[] A; delete[] b; delete[] y; delete[] cr;
		}
		#pragma endregion
	}

	for (size_t ispoke = 0; ispoke < num_spokes; ispoke++) delete[] spoke_end_points[ispoke];
	delete[] spoke_end_points; delete[] spoke_neighbors;  delete[] spoke_weights;
	delete[] neighbors; delete[] neighbor_hits;
	delete[] ff;

	return 0;
	#pragma endregion
}
