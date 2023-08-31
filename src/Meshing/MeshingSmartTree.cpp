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
//  MeshingSmartTree.cpp                                          Last modified (07/12/2022) //
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MeshingSmartTree.h"

const double DST_TOL = 1E-6;


MeshingSmartTree::MeshingSmartTree()
{
	init_global_variables();
	_xmin = new double[_num_dim];
	for (size_t idim = 0; idim < _num_dim; idim++) _xmin[idim] = DBL_MAX;
	_xmax = new double[_num_dim];
	for (size_t idim = 0; idim < _num_dim; idim++) _xmax[idim] = -DBL_MAX;
}

MeshingSmartTree::MeshingSmartTree(size_t num_dim)
{
	init_global_variables();
	_num_dim = num_dim;
	_xmin = new double[_num_dim];
	for (size_t idim = 0; idim < _num_dim; idim++) _xmin[idim] = DBL_MAX;
	_xmax = new double[_num_dim];
	for (size_t idim = 0; idim < _num_dim; idim++) _xmax[idim] = -DBL_MAX;
	
}

MeshingSmartTree::MeshingSmartTree(size_t num_dim, size_t num_points, double** points)
{
	init_global_variables();

	_num_dim = num_dim; _num_points = num_points; _points = points; _capacity = _num_points;

	_xmin = new double[_num_dim];
	for (size_t idim = 0; idim < _num_dim; idim++) _xmin[idim] = DBL_MAX;
	_xmax = new double[_num_dim];
	for (size_t idim = 0; idim < _num_dim; idim++) _xmax[idim] = -DBL_MAX;

	for (size_t ipoint = 0; ipoint < num_points; ipoint++)
	{
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			if (points[ipoint][idim] < _xmin[idim]) _xmin[idim] = points[ipoint][idim];
			if (points[ipoint][idim] > _xmax[idim]) _xmax[idim] = points[ipoint][idim];
		}
	}

	kd_tree_build_balanced();
}


MeshingSmartTree::~MeshingSmartTree()
{
	clear_memory();
}


int MeshingSmartTree::init_global_variables()
{
	/////////// kd-tree  Variables /////////////////////////////////////////////////
	_tree_left = 0; _tree_right = 0; _points = 0; _points_normal = 0; _points_attrib = 0; _marked = 0; _graph = 0;
	_num_points = 0; _marked_only = false; _num_points = 0; _capacity = 0;
	_num_dim = 3; _tree_origin = 0; _tree_max_height = 0; _auto_balance = true;

	_xmin = 0; _xmax = 0;

	reset_graph();
	return 0;
}

int MeshingSmartTree::init_random_samplers(size_t num_threads)
{
	#pragma region Init Random Samplers:
	_rsamplers.resize(num_threads);
	for (size_t thread_id = 0; thread_id < num_threads; thread_id++)
	{
		_rsamplers[thread_id] = new MeshingRandomSampler(thread_id);
	}
	return 0;
	#pragma endregion
}

int MeshingSmartTree::save_tree_csv(std::string file_name, size_t num_dim)
{
	#pragma region Save tree to CSV file:
	std::fstream file(file_name.c_str(), std::ios::out);
	// Spheres
	file << "x1coord";
	for (size_t idim = 1; idim < num_dim; idim++) file << ", x" << idim + 1 << "coord";

	for (size_t idim = 0; idim < _num_dim; idim++) file << ", n" << idim + 1;

	if (_points_attrib != 0 && _points_attrib[0] != 0)
	{
		size_t cap = _points_attrib[0][0];
		for (size_t j = 0; j < cap; j++) file << ", attrib" << j + 1;
	}
	file << std::endl;

	for (size_t i = 0; i < _num_points; i++)
	{
		file << std::setprecision(16) << _points[i][0];
		for (size_t idim = 1; idim < num_dim; idim++) file << ", " << _points[i][idim];

		if (_points_normal != 0 && _points_normal[i] != 0)
		{
			for (size_t idim = 0; idim < _num_dim; idim++) file << ", " << _points_normal[i][idim];
		}
		else
		{
			for (size_t idim = 0; idim < _num_dim; idim++) file << ", 0.0";
		}

		if (_points_attrib != 0 && _points_attrib[i] != 0)
		{
			size_t cap = _points_attrib[i][0];
			for (size_t j = 0; j < cap; j++) file << ", " << _points_attrib[i][j];
		}

		file << std::endl;
	}
	return 0;
	#pragma endregion
}

int MeshingSmartTree::re_enumerate_points_for_better_memory_access()
{
	size_t num_traversed(0);
	size_t* ordered_indices = new size_t[_num_points];
	kd_tree_get_nodes_order(0, _tree_origin, num_traversed, ordered_indices);

	size_t* point_new_index = new size_t[_num_points];
	for (size_t i = 0; i < _num_points; i++)
	{
		size_t current_point = ordered_indices[i];
		point_new_index[current_point] = i;
	}

	// new tree containers with the current order
	double** points = new double*[_capacity];
	double** points_normal = new double*[_capacity];
	size_t** points_attrib = new size_t*[_capacity];

	size_t* tree_left = new size_t[_capacity]; 
	size_t* tree_right = new size_t[_capacity];

	for (size_t i = 0; i < _num_points; i++)
	{
		size_t current_point = ordered_indices[i];
		points[i] = _points[current_point];
		points_normal[i] = _points_normal[current_point];
		points_attrib[i] = _points_attrib[current_point];

		tree_left[i] = point_new_index[_tree_left[current_point]];
		tree_right[i] = point_new_index[_tree_right[current_point]];
	}
	_tree_origin = point_new_index[_tree_origin];

	delete[] ordered_indices; delete[] point_new_index;
	delete[] _points; delete[] _points_normal; delete[] _points_attrib;
	delete[] _tree_left; delete[] _tree_right;

	_points = points; _points_normal = points_normal; _points_attrib = points_attrib;
	_tree_left = tree_left; _tree_right = tree_right;

	return 0;
}

int MeshingSmartTree::kd_tree_get_nodes_order(size_t d_index, size_t node_index,                                // indices to traverse the kd-tree
	                                   size_t& num_traversed, size_t* ordered_indices)
{
	ordered_indices[num_traversed] = node_index; num_traversed++;
	if (_tree_right[node_index] != node_index) kd_tree_get_nodes_order(d_index + 1, _tree_right[node_index], num_traversed, ordered_indices);
	if (_tree_left[node_index] != node_index) kd_tree_get_nodes_order(d_index + 1, _tree_left[node_index], num_traversed, ordered_indices);
	return 0;
}

int MeshingSmartTree::get_tokens(std::string line, char separator, std::vector<std::string>& tokens)
{
	tokens.clear();
	size_t start(0);
	while (true)
	{
		size_t end = line.find_first_of(separator, start);
		if (end == std::string::npos) end = line.size();

		tokens.push_back(line.substr(start, end - start));
		start = end + 1;
		
		if (end == line.size())
		{
			// end of string, exit loop
			break;
		}
	}
	return 0;
}

int MeshingSmartTree::load_tree_csv(std::string file_name, size_t num_dim)
{
	#pragma region Load tree from CSV file:
	clock_t start_time, end_time; double cpu_time;
	start_time = clock();

	// vcm_cout << "VoroCrust::Loading a Smart kd-tree from csv file " << file_name.c_str() << ":" << std::endl;

	//open file
	std::ifstream tmpfile(file_name.c_str());

	if (tmpfile.is_open() && tmpfile.good())
	{
		std::string line = "";
		size_t iline(0), num_entries(0), num_attrib(0);
		double* x(0); double* normal(0); size_t* attrib(0);
		while (getline(tmpfile, line))
		{
			std::istringstream iss(line);
			std::vector<std::string> tokens;
			get_tokens(line, ',', tokens);

			if (iline == 0)
			{
				num_entries = tokens.size();
				x = new double[num_dim];
				if (num_entries >= num_dim + _num_dim) normal = new double[_num_dim];
				if (num_entries > num_dim + _num_dim)
				{
					num_attrib = num_entries - num_dim - _num_dim;
					attrib = new size_t[num_attrib];
				}
				iline++;
				continue;
			}
			if (tokens.size() == 0) continue;

			for (size_t i = 0; i < num_dim; i++) x[i] = _memo.string_to_double(tokens[i]);

			if (num_entries >= num_dim + _num_dim)
			{
				for (size_t i = 0; i < _num_dim; i++) normal[i] = _memo.string_to_double(tokens[num_dim + i]);
			}

			if (num_entries > num_dim + _num_dim)
			{
				for (size_t i = 0; i < num_attrib; i++) attrib[i] = _memo.string_to_size_t(tokens[num_dim +_num_dim + i]);
			}
			add_tree_point(num_dim, x, normal, attrib);
			iline++;
		}
		delete[] x;
		if (normal != 0) delete[] normal;
		if (attrib != 0) delete[] attrib;
	}
	end_time = clock();
	cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
	// vcm_cout << "  * executed in " << cpu_time << " seconds!" << std::endl;
	// vcm_cout << "  * Number of tree points = " << _num_points << std::endl;
	return 0;
	#pragma endregion
}

int MeshingSmartTree::reset_graph()
{
	#pragma region Reset Point Graph:

	if (_num_points == 0) return 0;

	if (_graph != 0)
	{
		for (size_t i = 0; i < _graph_cap; i++)
		{
			if (_graph[i] != 0) delete[] _graph[i];
		}
		delete[] _graph;
	}
	_graph = 0;
	return 0;
	#pragma endregion
}

int MeshingSmartTree::graph_connect_nodes(size_t ipoint, size_t jpoint)
{
	graph_connect(ipoint, jpoint);
	graph_connect(jpoint, ipoint);
	return 0;
}

int MeshingSmartTree::graph_connect(size_t ipoint, size_t jpoint)
{
	#pragma region Add Graph Connection:

	if (ipoint >= _num_points) return 1;

	if (_graph == 0)
	{
		_graph_cap = _capacity;
		_graph = new size_t*[_graph_cap];
		for (size_t i = 0; i < _graph_cap; i++) _graph[i] = 0;
	}
	else if (_graph_cap < _capacity)
	{
		// increase the size of the graph list
		size_t** graph_new = new size_t*[_capacity];
		for (size_t i = 0; i < _graph_cap; i++) graph_new[i] = _graph[i];
		for (size_t i = _graph_cap; i < _capacity; i++) graph_new[i] = 0;
		delete[] _graph;
		_graph = graph_new;
		_graph_cap = _capacity;
	}

	if (_graph[ipoint] == 0)
	{
		_graph[ipoint] = new size_t[12];
		_graph[ipoint][0] = 10; _graph[ipoint][1] = 0;
	}

	if (graph_connected(ipoint, jpoint)) return 1;

	size_t edges_cap(_graph[ipoint][0]), num_edges(_graph[ipoint][1]);

	_graph[ipoint][2 + num_edges] = jpoint;
	_graph[ipoint][1]++; num_edges++;

	if (num_edges == edges_cap)
	{
		size_t edges_cap_new = 2 * edges_cap;
		size_t* edges_new = new size_t[edges_cap_new + 2];
		for (size_t i = 0; i < num_edges; i++) edges_new[2 + i] = _graph[ipoint][2 + i];
		edges_new[0] = edges_cap_new; edges_new[1] = num_edges;
		delete[] _graph[ipoint];
		_graph[ipoint] = edges_new;
	}
	return 0;
	#pragma endregion
}

bool MeshingSmartTree::graph_connected(size_t ipoint, size_t jpoint)
{
	#pragma region Connect Graph points:
	if (_graph == 0)          return false;
	if (ipoint >= _graph_cap) return false;
	if (_graph[ipoint] == 0)  return false;

	size_t edges_cap(_graph[ipoint][0]), num_edges(_graph[ipoint][1]);
	for (size_t ied = 0; ied < num_edges; ied++)
	{
		if (_graph[ipoint][2 + ied] == jpoint) return true; // connection already exist
	}
	return false;
	#pragma endregion
}

int MeshingSmartTree::add_tree_point(size_t num_dim, double* x, double* normal, size_t* attrib)
{
	#pragma region Add a tree point:

	if (_xmin == 0)
	{
		_xmin = new double[_num_dim];
		for (size_t idim = 0; idim < _num_dim; idim++) _xmin[idim] = DBL_MAX;
	}

	if (_xmax == 0)
	{
		_xmax = new double[_num_dim];
		for (size_t idim = 0; idim < _num_dim; idim++) _xmax[idim] = -DBL_MAX;
	}

	for (size_t idim = 0; idim < _num_dim; idim++)
	{
		_xmin[idim] = fmin(_xmin[idim], x[idim]);
		_xmax[idim] = fmax(_xmax[idim], x[idim]);
	}

	if (_num_points == 0)
	{
		_capacity = 100;
		_points = new double*[_capacity];
		_points_normal = new double*[_capacity];
		_points[0] = new double[num_dim];
		_points_normal[0] = new double[num_dim];
		for (size_t idim = 0; idim < num_dim; idim++)
		{
			_points[0][idim] = x[idim];
			if (normal == 0)                      _points_normal[0][idim] = 0.0;
			else if (idim < _num_dim)             _points_normal[0][idim] = normal[idim];
		}

		_points_attrib = new size_t*[_capacity];
		if (attrib == 0) _points_attrib[0] = 0;
		else
		{
			size_t num_attrib = attrib[0];
			_points_attrib[0] = new size_t[num_attrib];
			for (size_t i = 0; i < num_attrib; i++) _points_attrib[0][i] = attrib[i];
		}

		_tree_left = new size_t[_capacity]; _tree_right = new size_t[_capacity];
		_tree_origin = 0; _tree_left[0] = 0; _tree_right[0] = 0;
		_tree_max_height = 1;
		_num_points++;
	}
	else
	{
		_points[_num_points] = new double[num_dim];
		_points_normal[_num_points] = new double[num_dim];

		for (size_t idim = 0; idim < num_dim; idim++) _points[_num_points][idim] = x[idim];

		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			if (normal == 0) _points_normal[_num_points][idim] = 0.0;
			else             _points_normal[_num_points][idim] = normal[idim];
		}

		if (attrib == 0) _points_attrib[_num_points] = 0;
		else
		{
			size_t num_attrib = attrib[0];
			_points_attrib[_num_points] = new size_t[num_attrib];
			for (size_t i = 0; i < num_attrib; i++) _points_attrib[_num_points][i] = attrib[i];
		}

		_tree_left[_num_points] = _num_points; _tree_right[_num_points] = _num_points;

		kd_tree_add_point(_num_points);

		_num_points++;

		if (_num_points == _capacity)
		{
			size_t capacity = _capacity * 2;
			double** points = new double*[capacity];
			double** points_normal = new double*[capacity];
			size_t** points_attrib = new size_t*[capacity];
			size_t* tree_left = new size_t[capacity];
			size_t* tree_right = new size_t[capacity];
			for (size_t i = 0; i < _capacity; i++)
			{
				points[i] = _points[i];
				points_normal[i] = _points_normal[i];
				points_attrib[i] = _points_attrib[i];
				tree_left[i] = _tree_left[i];
				tree_right[i] = _tree_right[i];
			}
			delete[] _points;        _points = points;
			delete[] _points_normal; _points_normal = points_normal;
			delete[] _points_attrib; _points_attrib = points_attrib;
			delete[] _tree_left;     _tree_left = tree_left;
			delete[] _tree_right;    _tree_right = tree_right;
			_capacity = capacity;
		}

		if (_auto_balance &&  _tree_max_height > log2(1 + _num_points) * 10)
		{
			kd_tree_build_balanced();
		}
	}
	return 0;
	#pragma endregion
}

int MeshingSmartTree::get_tree_point(size_t point_index, double* x)
{
	for (size_t idim = 0; idim < _num_dim; idim++) x[idim] = _points[point_index][idim];
	return 0;
}

int MeshingSmartTree::get_tree_point(size_t point_index, size_t num_dim, double* x)
{
	for (size_t idim = 0; idim < num_dim; idim++) x[idim] = _points[point_index][idim];
	return 0;
}

bool MeshingSmartTree::get_tree_point_attrib(size_t point_index, size_t attrib_index, size_t &point_attrib)
{
	if (_points_attrib[point_index] == 0) return false;
	if (attrib_index >= _points_attrib[point_index][0]) return false;
	point_attrib = _points_attrib[point_index][1 + attrib_index];
	return true;
}

double MeshingSmartTree::get_tree_point_attrib(size_t point_index, size_t attrib_index)
{
	return _points[point_index][_num_dim + attrib_index];
}

int MeshingSmartTree::set_tree_point_attrib(size_t point_index, size_t attrib_index, size_t attrib)
{
	_points_attrib[point_index][1 + attrib_index] = attrib;
	return 0;
}

int MeshingSmartTree::set_tree_point_attrib(size_t point_index, size_t attrib_index, double attrib)
{
	_points[point_index][_num_dim + attrib_index] = attrib;
	return 0;
}

void MeshingSmartTree::clear_memory()
{
	#pragma region Clear Memory:
	if (_tree_left != 0)
	{
		delete[] _tree_left; _tree_left = 0;
	}
	if (_tree_right != 0)
	{
		delete[] _tree_right; _tree_right = 0;
	}
	if (_points_attrib != 0)
	{
		for (size_t ipoint = 0; ipoint < _num_points; ipoint++) delete[] _points_attrib[ipoint];
		delete[] _points_attrib; _points_attrib = 0;
	}
	if (_points != 0)
	{
		for (size_t ipoint = 0; ipoint < _num_points; ipoint++) delete[] _points[ipoint];
		delete[] _points; _points = 0;
	}
	if (_points_normal != 0)
	{
		for (size_t ipoint = 0; ipoint < _num_points; ipoint++) delete[] _points_normal[ipoint];
		delete[] _points_normal; _points_normal = 0;
	}

	if (_marked != 0)
	{
		delete[] _marked;
		_marked = 0;
		_marked_only = false;
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

	reset_graph();


	size_t num_threads(_rsamplers.size());
	for (size_t thread_id = 0; thread_id < num_threads; thread_id++)
	{
		delete _rsamplers[thread_id];
	}

	_tree_max_height = 0;
	_num_points = 0; _capacity = 0;
	#pragma endregion
}

int MeshingSmartTree::kd_tree_build_balanced()
{
	#pragma region Build Balanced kd-tree:
	//vcm_cout << "SmartkdTree::Building_Balanced_Tree:" << std::endl;

	clock_t start_time, end_time; double cpu_time;
	start_time = clock();

	if (_tree_left == 0) _tree_left = new size_t[_capacity];
	if (_tree_right == 0) _tree_right = new size_t[_capacity];

	size_t* tree_nodes_sorted = new size_t[_num_points];
	for (size_t i = 0; i < _num_points; i++) tree_nodes_sorted[i] = i;

	for (size_t iseed = 0; iseed < _num_points; iseed++)
	{
		_tree_left[iseed] = iseed; _tree_right[iseed] = iseed;
	}
	_tree_origin = _capacity;

	size_t target_pos = _num_points / 2; _tree_max_height = 0;
	kd_tree_balance_quicksort(target_pos, 0, _num_points - 1, 0, tree_nodes_sorted);

	delete[] tree_nodes_sorted;

	end_time = clock();
	cpu_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

	//vcm_cout << "      * Executed in " << cpu_time << " seconds." << std::endl;
	//vcm_cout << "      * Number of tree levels = " << _tree_max_height << std::endl;
	//vcm_cout << "      * Number of perfectly balanced tree levels = " << ceil(log2(1 + _num_points)) << std::endl;

	return 0;
	#pragma endregion
}

int MeshingSmartTree::get_closest_tree_point(size_t tree_point_index, size_t &closest_tree_point, double &closest_distance)
{
	#pragma region Closest Neighbor Search using tree:

	closest_tree_point = _num_points;
	kd_tree_get_closest_seed(tree_point_index, 0, _tree_origin, closest_tree_point, closest_distance);

	return 0;
	#pragma endregion
}

int MeshingSmartTree::get_closest_tree_point(double* x, size_t &closest_tree_point, double &closest_distance)
{
	#pragma region Closest Neighbor Search using tree:
	closest_tree_point = _num_points;
	size_t num_nodes_visited = 0;
	if (_num_points == 0) return 1;
	kd_tree_get_closest_seed(x, 0, _tree_origin, closest_tree_point, closest_distance, num_nodes_visited);
	return 0;
	#pragma endregion
}

int MeshingSmartTree::get_closest_tree_point(double* x, size_t num_exculded_tree_points, size_t* exculded_tree_points, size_t &closest_tree_point, double &closest_distance)
{
	#pragma region Closest Neighbor Search using a tree:

	closest_tree_point = _num_points;
	kd_tree_get_closest_seed(x, 0, _tree_origin, num_exculded_tree_points, exculded_tree_points, closest_tree_point, closest_distance);
	return 0;

	#pragma endregion
}

// used for estimating sizing function in VC_Wild
int MeshingSmartTree::get_closest_non_smooth_tree_point(double* x, double feature_TOL, double smoothness_threshold, size_t &closest_tree_point, double &closest_distance)
{
	closest_tree_point = _num_points;
	kd_tree_get_closest_non_smooth_seed(x, feature_TOL, smoothness_threshold, 0, _tree_origin, closest_tree_point, closest_distance);
	return 0;
}

// used for estimating sizing function in VC_Wild
int MeshingSmartTree::get_closest_non_smooth_tree_point(double* x, double* x_normal, double smoothness_threshold, size_t& closest_tree_point, double& closest_distance)
{
	closest_tree_point = _num_points;
	kd_tree_get_closest_non_smooth_seed(x, x_normal, smoothness_threshold, 0, _tree_origin, closest_tree_point, closest_distance);
	return 0;
}


// used for detecting sharp corners/edges in VC_Wild
int MeshingSmartTree::get_closest_non_smooth_tree_point(double* x, size_t num_normals, double** normals, double smoothness_threshold, size_t &closest_tree_point, double &closest_distance)
{
	closest_tree_point = _num_points;
	kd_tree_get_closest_non_smooth_seed(x, num_normals, normals, smoothness_threshold, 0, _tree_origin, closest_tree_point, closest_distance);
	return 0;
}


// used for estimating curves sizing function in VC_Wild
int MeshingSmartTree::get_closest_non_smooth_tree_edge_point(double* x, double feature_TOL, double smoothness_threshold, size_t &closest_tree_point, double &closest_distance)
{
	closest_tree_point = _num_points;
	kd_tree_get_closest_non_smooth_edge_seed(x, feature_TOL, smoothness_threshold, 0, _tree_origin, closest_tree_point, closest_distance);
	return 0;
}

int MeshingSmartTree::get_tree_points_in_sphere(double* x, double r, size_t &num_points_in_sphere, size_t* &points_in_sphere, size_t &capacity)
{
	#pragma region tree neighbor search:
	num_points_in_sphere = 0;
	kd_tree_get_seeds_in_sphere(x, r, 0, _tree_origin, num_points_in_sphere, points_in_sphere, capacity);
	return 0;
	#pragma endregion
}

int MeshingSmartTree::get_Voronoi_neighbors(size_t num_spokes, size_t num_layers, size_t num_dim, double* x,
	                                 size_t &num_neighbors, size_t* &neighbors, size_t* &neighbor_hits, size_t* &neighbor_level)
{
	get_Voronoi_neighbors(num_spokes, num_layers, num_dim, x, 0, 0, num_neighbors, neighbors, neighbor_hits, neighbor_level);
	return 0;
}

int MeshingSmartTree::get_Voronoi_neighbors(size_t num_spokes, size_t num_layers, size_t num_dim, double* x,
									 double** spoke_points, size_t* spoke_neighbors,
	                                 size_t &num_neighbors, size_t* &neighbors, size_t* &neighbor_hits)
{
	size_t* neighbor_level(0);
	get_Voronoi_neighbors(num_spokes, num_layers, num_dim, x, spoke_points, spoke_neighbors, num_neighbors, neighbors, neighbor_hits, neighbor_level);
	if (neighbor_level != 0) delete[] neighbor_level;
	return 0;
}

int MeshingSmartTree::get_Voronoi_neighbors(size_t num_desired_neighbors, size_t num_dim, double* x,
	                                 size_t &num_neighbors, size_t* &neighbors)
{
    #pragma region Get Voronoi Neighbors:

	num_neighbors = 0;

	size_t cap(100);
	neighbors = new size_t[cap];

	size_t iclosest; double hclosest(DBL_MAX);
	get_closest_tree_point(x, 0, 0, iclosest, hclosest);

	if (hclosest < 1E-10)
	{
		// query point is a seed
		neighbors[0] = iclosest;
		num_neighbors = 1;

		hclosest = DBL_MAX;
		get_closest_tree_point(x, num_neighbors, neighbors, iclosest, hclosest);
	}

	double* spoke_dir = new double[num_dim];
	while (num_neighbors< num_desired_neighbors)
	{
		_rsampler.sample_uniformly_from_unit_sphere(spoke_dir, num_dim);

		double spoke_length(DBL_MAX); size_t neighbor_seed(_num_points);
		trim_spoke(x, num_neighbors, neighbors, x, spoke_dir, spoke_length, neighbor_seed);

		if (spoke_length == DBL_MAX) continue; // An unbounded spoke

		neighbors[num_neighbors] = neighbor_seed;
		num_neighbors++;
		if (num_neighbors == cap)
		{
			#pragma region Expand neighbors containers:
			size_t new_cap(cap * 2);
			size_t* new_neighbors = new size_t[new_cap];
			for (size_t i = 0; i < cap; i++)
			{
				new_neighbors[i] = neighbors[i];
			}
			cap = new_cap;
			delete[] neighbors;
			neighbors = new_neighbors;
			#pragma endregion
		}
	}
	delete[] spoke_dir;

	return 0;
	#pragma endregion
}

int MeshingSmartTree::get_Voronoi_neighbors(size_t num_spokes, size_t num_layers, size_t num_dim, double* x,
	                                 double** spoke_points, size_t* spoke_neighbors,
	                                 size_t &num_neighbors, size_t* &neighbors, size_t* &neighbor_hits, size_t* &neighbor_level)
{
	#pragma region Get Voronoi Neighbors:

	num_neighbors = 0;

	size_t cap(100);
	neighbors = new size_t[cap];
	neighbor_hits = new size_t[cap];
	neighbor_level = new size_t[cap];

	size_t iclosest; double hclosest(DBL_MAX);
	get_closest_tree_point(x, 0, 0, iclosest, hclosest);

	size_t num_old_neighbors(0);
	size_t* old_neighbors(0);

	if (hclosest < 1E-10)
	{
		// layer -1: query point is a seed
		neighbors[0] = iclosest; neighbor_hits[0] = SIZE_MAX; neighbor_level[0] = 0;
		num_neighbors = 1;

		num_old_neighbors = 1;
		old_neighbors = new size_t[1];
		old_neighbors[0] = iclosest;
		hclosest = DBL_MAX;
		get_closest_tree_point(x, num_old_neighbors, old_neighbors, iclosest, hclosest);
	}


	for (size_t ilayer = 0; ilayer < num_layers; ilayer++)
	{
		for (size_t ispoke = 0; ispoke < num_spokes; ispoke++)
		{
			double* spoke_dir = new double[num_dim];
			_rsampler.sample_uniformly_from_unit_sphere(spoke_dir, num_dim);

			double spoke_length(DBL_MAX); size_t neighbor_seed(_num_points);
			trim_spoke(x, num_old_neighbors, old_neighbors, x, spoke_dir, spoke_length, neighbor_seed);

			if (spoke_points != 0)
			{
				spoke_points[ispoke] = spoke_dir;
				spoke_neighbors[ispoke] = SIZE_MAX;
			}
			else delete[] spoke_dir;

			if (spoke_length == DBL_MAX) continue; // An unbounded spoke

			if (spoke_points != 0)
			{
				for (size_t idim = 0; idim < num_dim; idim++) spoke_points[ispoke][idim] = x[idim] + spoke_length * spoke_dir[idim];
				spoke_neighbors[ispoke] = neighbor_seed;
			}

			bool found(false);
			for (size_t i = 0; i < num_neighbors; i++)
			{
				if (neighbors[i] == neighbor_seed)
				{
					neighbor_hits[i]++;
					found = true;
					break;
				}
			}

			if (!found)
			{
				neighbors[num_neighbors] = neighbor_seed;
				neighbor_hits[num_neighbors] = 1;
				neighbor_level[num_neighbors] = ilayer + 1;
				num_neighbors++;
				if (num_neighbors == cap)
				{
					#pragma region Expand neighbors containers:
					size_t new_cap(cap * 2);
					size_t* new_neighbors = new size_t[new_cap];
					size_t* new_neighbors_hit = new size_t[new_cap];
					size_t* new_neighbor_level = new size_t[new_cap];
					for (size_t i = 0; i < cap; i++)
					{
						new_neighbors[i] = neighbors[i];
						new_neighbors_hit[i] = neighbor_hits[i];
						new_neighbor_level[i] = neighbor_level[i];
					}
					cap = new_cap;
					delete[] neighbors;
					delete[] neighbor_hits;
					delete[] neighbor_level;
					neighbors = new_neighbors;
					neighbor_hits = new_neighbors_hit;
					neighbor_level = new_neighbor_level;
					#pragma endregion
				}
			}
		}

		if (old_neighbors != 0)
		{
			delete[] old_neighbors;
			old_neighbors = 0;
		}

		if (ilayer < num_layers - 1)
		{
			num_old_neighbors = num_neighbors;
			old_neighbors = new size_t[num_neighbors];
			for (size_t ii = 0; ii < num_neighbors; ii++) old_neighbors[ii] = neighbors[ii];
		}
	}

	if (old_neighbors != 0) delete[] old_neighbors;

	return 0;
	#pragma endregion
}


int MeshingSmartTree::get_Voronoi_neighbors(size_t thread_id, size_t num_spokes, size_t num_layers, size_t num_dim, double* x,
	                                 double** spoke_points, size_t* spoke_neighbors,
	                                 size_t& num_neighbors, size_t*& neighbors, size_t*& neighbor_hits)
{
	size_t* neighbor_level(0);
	get_Voronoi_neighbors(thread_id, num_spokes, num_layers, num_dim, x, spoke_points, spoke_neighbors, num_neighbors, neighbors, neighbor_hits, neighbor_level);
	if (neighbor_level != 0) delete[] neighbor_level;
	return 0;
}

int MeshingSmartTree::get_Voronoi_neighbors(size_t thread_id, size_t num_spokes, size_t num_layers, size_t num_dim, double* x,
	                                 double** spoke_points, size_t* spoke_neighbors,
	                                 size_t& num_neighbors, size_t*& neighbors, size_t*& neighbor_hits, size_t*& neighbor_level)
{
	#pragma region Get Voronoi Neighbors:

	num_neighbors = 0;

	size_t cap(100);
	neighbors = new size_t[cap];
	neighbor_hits = new size_t[cap];
	neighbor_level = new size_t[cap];

	size_t iclosest; double hclosest(DBL_MAX);
	get_closest_tree_point(x, 0, 0, iclosest, hclosest);

	size_t num_old_neighbors(0);
	size_t* old_neighbors(0);

	if (hclosest < 1E-10)
	{
		// layer -1: query point is a seed
		neighbors[0] = iclosest; neighbor_hits[0] = SIZE_MAX; neighbor_level[0] = 0;
		num_neighbors = 1;

		num_old_neighbors = 1;
		old_neighbors = new size_t[1];
		old_neighbors[0] = iclosest;
		hclosest = DBL_MAX;
		get_closest_tree_point(x, num_old_neighbors, old_neighbors, iclosest, hclosest);
	}


	for (size_t ilayer = 0; ilayer < num_layers; ilayer++)
	{
		for (size_t ispoke = 0; ispoke < num_spokes; ispoke++)
		{
			double* spoke_dir = new double[num_dim];
			_rsamplers[thread_id]->sample_uniformly_from_unit_sphere(spoke_dir, num_dim);

			double spoke_length(DBL_MAX); size_t neighbor_seed(_num_points);
			trim_spoke(x, num_old_neighbors, old_neighbors, x, spoke_dir, spoke_length, neighbor_seed);

			if (spoke_points != 0)
			{
				spoke_points[ispoke] = spoke_dir;
				spoke_neighbors[ispoke] = SIZE_MAX;
			}
			else delete[] spoke_dir;

			if (spoke_length == DBL_MAX) continue; // Anunbounded spoke

			if (spoke_points != 0)
			{
				for (size_t idim = 0; idim < num_dim; idim++) spoke_points[ispoke][idim] = x[idim] + spoke_length * spoke_dir[idim];
				spoke_neighbors[ispoke] = neighbor_seed;
			}

			bool found(false);
			for (size_t i = 0; i < num_neighbors; i++)
			{
				if (neighbors[i] == neighbor_seed)
				{
					neighbor_hits[i]++;
					found = true;
					break;
				}
			}

			if (!found)
			{
				neighbors[num_neighbors] = neighbor_seed;
				neighbor_hits[num_neighbors] = 1;
				neighbor_level[num_neighbors] = ilayer + 1;
				num_neighbors++;
				if (num_neighbors == cap)
				{
					#pragma region Expand neighbors containers:
					size_t new_cap(cap * 2);
					size_t* new_neighbors = new size_t[new_cap];
					size_t* new_neighbors_hit = new size_t[new_cap];
					size_t* new_neighbor_level = new size_t[new_cap];
					for (size_t i = 0; i < cap; i++)
					{
						new_neighbors[i] = neighbors[i];
						new_neighbors_hit[i] = neighbor_hits[i];
						new_neighbor_level[i] = neighbor_level[i];
					}
					cap = new_cap;
					delete[] neighbors;
					delete[] neighbor_hits;
					delete[] neighbor_level;
					neighbors = new_neighbors;
					neighbor_hits = new_neighbors_hit;
					neighbor_level = new_neighbor_level;
					#pragma endregion
				}
			}
		}

		if (old_neighbors != 0)
		{
			delete[] old_neighbors;
			old_neighbors = 0;
		}

		if (ilayer < num_layers - 1)
		{
			num_old_neighbors = num_neighbors;
			old_neighbors = new size_t[num_neighbors];
			for (size_t ii = 0; ii < num_neighbors; ii++) old_neighbors[ii] = neighbors[ii];
		}
	}

	if (old_neighbors != 0) delete[] old_neighbors;

	return 0;
	#pragma endregion
}

int MeshingSmartTree::trim_spoke(double* cell_seed, size_t num_involved_seeds, size_t* involved_seeds, double* spoke_origin, double* spoke_dir, double &spoke_length, size_t &neighbor_seed)
{
	#pragma region Trim a spoke using Tree points:

	if (num_involved_seeds == 0)
	{
		//vcm_cout << "Error:: Spoke trimming require ate least one seed to be involved!" << std::endl;
		//return 1;
	}

	double origin_dst(0.0);
	for (size_t idim = 0; idim < _num_dim; idim++)
	{
		double dx = spoke_origin[idim] - cell_seed[idim];
		origin_dst += dx * dx;
	}
	origin_dst = sqrt(origin_dst);

	kd_tree_trim_spoke(cell_seed, num_involved_seeds, involved_seeds, spoke_origin, spoke_dir, origin_dst, 0, _tree_origin, spoke_length, neighbor_seed);

	if (num_involved_seeds == 0) return 0;


	// verify the result
	double* xend = new double[_num_dim];
	for (size_t idim = 0; idim < _num_dim; idim++) xend[idim] = spoke_origin[idim] + spoke_length * spoke_dir[idim];

	double dst_seed(0.0);
	for (size_t idim = 0; idim < _num_dim; idim++)
	{
		double dx = xend[idim] - cell_seed[idim];
		dst_seed += dx * dx;
	}
	dst_seed = sqrt(dst_seed);

	size_t iclosest; double hclosest(DBL_MAX);
	get_closest_tree_point(xend, iclosest, hclosest);

	bool found(false);
	for (size_t i = 0; i < num_involved_seeds; i++)
	{
		if (iclosest == involved_seeds[i])
		{
			found = true;
			break;
		}
	}

	if (!found && hclosest < dst_seed - 1E-6)
	{
		size_t jclosest; double hhclosest(DBL_MAX);
		get_closest_tree_point(spoke_origin, jclosest, hhclosest);

		spoke_length = DBL_MAX;

		kd_tree_trim_spoke(cell_seed, num_involved_seeds, involved_seeds, spoke_origin, spoke_dir, origin_dst, 0, _tree_origin, spoke_length, neighbor_seed);
	}
	delete[] xend;

	return 0;
	#pragma endregion
}

int MeshingSmartTree::get_3d_Voronoi_cell(size_t seed_index, size_t &num_cell_neighbors, size_t* &cell_neighbors, double &cell_volume,
	                               size_t* &num_faces_corners, double*** &faces_corners, double* &faces_area)
{
	#pragma region Construct 3d Voronoi Cell:
	if (_num_dim != 3) return 1;

	size_t num_spokes(1);

	cell_volume = 0.0;

	num_cell_neighbors = 0;
	size_t cell_neighbors_cap(100);
	cell_neighbors = new size_t[cell_neighbors_cap];
	num_faces_corners = new size_t[cell_neighbors_cap];
	faces_corners = new double**[cell_neighbors_cap];
	faces_area = new double[cell_neighbors_cap];

	size_t num_additional_neighbors(0), additional_neighbors_cap(100);
	size_t* additional_neighbors = new size_t[additional_neighbors_cap];
	double** additional_neighbors_x = new double*[additional_neighbors_cap];

	double* spoke_dir = new double[3];
	double* xface = new double[3];

	size_t* seeds = new size_t[2]; seeds[0] = seed_index;

	size_t ispoke(0);
	size_t closest_seed;  double hclosest(DBL_MAX);
	bool unbounded_cell(false);
	while (true)
	{
		ispoke++;
		size_t neighbor_seed(seed_index);
		if (ispoke == 1)
		{
			#pragma region Voronoi face with closest neighbor:
			get_closest_tree_point(seed_index, closest_seed, hclosest);
			for (size_t idim = 0; idim < 3; idim++) xface[idim] = 0.5 * (_points[seed_index][idim] + _points[closest_seed][idim]);
			neighbor_seed = closest_seed;
			#pragma endregion
		}
		else if (ispoke <= num_spokes)
		{
			#pragma region A random Voronoi Face using VoroSpoke:
			_rsampler.sample_uniformly_from_unit_sphere(spoke_dir, 3);

			double spoke_length(DBL_MAX);
			trim_spoke(_points[seed_index], 1, seeds, _points[seed_index], spoke_dir, spoke_length, neighbor_seed);

			if (spoke_length == DBL_MAX)
			{
				// an infinit direction
				unbounded_cell = true;
				continue;
			}
			for (size_t idim = 0; idim < _num_dim; idim++) xface[idim] = _points[seed_index][idim] + spoke_length * spoke_dir[idim];

			bool prior_neighbor(false);
			for (size_t i = 0; i < num_cell_neighbors; i++)
			{
				if (cell_neighbors[i] == neighbor_seed)
				{
					prior_neighbor = true;
					break;
				}
			}
			if (prior_neighbor) continue; // We already constructed the associated face
			#pragma endregion
		}
		else if (num_additional_neighbors > 0)
		{
			#pragma region remaining Voronoi Faces:
			size_t neighbor_index = num_additional_neighbors - 1;
			neighbor_seed = additional_neighbors[neighbor_index];
			for (size_t idim = 0; idim < 3; idim++)  xface[idim] = additional_neighbors_x[neighbor_index][idim];
			num_additional_neighbors--;
			#pragma endregion
		}
		else break;

		// Construct a Voronoi Face between seed_index and neighbor_seed
		seeds[1] = neighbor_seed;
		double* si = _points[seed_index]; double* sj = _points[neighbor_seed];

		size_t num_face_corners(0);
		double** face_corners(0);
		double face_area(0.0), face_height(0.0);

		get_3d_Voronoi_face(seed_index, neighbor_seed, xface, num_face_corners, face_corners, face_area, face_height);

		for (size_t icorner = 0; icorner < num_face_corners; icorner++)
		{
			#pragma region Collecting additional Neighbors:
			double r(0.0);
			for (size_t idim = 0; idim < 3; idim++)
			{
				double dx = face_corners[icorner][idim] - _points[seed_index][idim];
				r += dx * dx;
			}
			r = sqrt(r);

			size_t num_neighbors(0), neighbors_cap(100);
			size_t* neighbors = new size_t[neighbors_cap];
			get_tree_points_in_sphere(face_corners[icorner], r, num_neighbors, neighbors, neighbors_cap);
			for (size_t ii = 0; ii < num_neighbors; ii++)
			{
				size_t neighbor_index = neighbors[ii];
				if (neighbor_index == seed_index) continue;
				if (neighbor_index == neighbor_seed) continue;

				bool found(false);
				for (size_t j = 0; j < num_cell_neighbors; j++)
				{
					if (neighbor_index == cell_neighbors[j])
					{
						found = true; break;
					}
				}
				if (found) continue; // A prior neighbor

				for (size_t j = 0; j < num_additional_neighbors; j++)
				{
					if (neighbor_index == additional_neighbors[j])
					{
						found = true; break;
					}
				}
				if (found) continue; // already exists in the additional neighbors list

				additional_neighbors[num_additional_neighbors] = neighbor_index;
				additional_neighbors_x[num_additional_neighbors] = face_corners[icorner];
				num_additional_neighbors++;

				if (num_additional_neighbors == additional_neighbors_cap)
				{
					#pragma region Increase Containers capacity:
					size_t new_cap = 2 * additional_neighbors_cap;
					size_t* new_additional_neighbors = new size_t[new_cap];
					double** new_additional_neighbors_x = new double* [new_cap];
					for (size_t i = 0; i < num_additional_neighbors; i++)
					{
						new_additional_neighbors[i] = additional_neighbors[i];
						new_additional_neighbors_x[i] = additional_neighbors_x[i];
					}

					delete[] additional_neighbors; delete[] additional_neighbors_x;
					additional_neighbors = new_additional_neighbors;
					additional_neighbors_x = new_additional_neighbors_x;
					additional_neighbors_cap = new_cap;
					#pragma endregion
				}
			}
			delete[] neighbors;
			#pragma endregion
		}

		cell_neighbors[num_cell_neighbors] = neighbor_seed;
		for (size_t i = 0; i < num_additional_neighbors; i++)
		{
			if (additional_neighbors[i] == neighbor_seed)
			{
				// remove neighbor from the additional neighbors list
				additional_neighbors[i] = additional_neighbors[num_additional_neighbors - 1];
				additional_neighbors_x[i] = additional_neighbors_x[num_additional_neighbors - 1];
				num_additional_neighbors--;
			}
		}

		num_faces_corners[num_cell_neighbors] = num_face_corners;
		faces_corners[num_cell_neighbors] = face_corners;
		faces_area[num_cell_neighbors] = face_area;

		if (face_area == DBL_MAX)
		{
			// unbounded face
			unbounded_cell = true; cell_volume = DBL_MAX;
		}
		else
		{
			// a non-degenerate face
			cell_volume += faces_area[num_cell_neighbors] * face_height / 3;
		}

		num_cell_neighbors++;

		if (num_cell_neighbors == cell_neighbors_cap)
		{
			#pragma region Increase Containers capacity:
			size_t new_cap = 2 * cell_neighbors_cap;
			size_t* new_cell_neighbors = new size_t[new_cap];
			size_t* new_num_faces_corners = new size_t[new_cap];
			double*** new_faces_corners = new double** [new_cap];
			double* new_faces_area = new double[new_cap];
			for (size_t i = 0; i < num_cell_neighbors; i++)
			{
				new_cell_neighbors[i] = cell_neighbors[i];
				new_num_faces_corners[i] = num_faces_corners[i];
				new_faces_corners[i] = faces_corners[i];
				new_faces_area[i] = faces_area[i];
			}
			delete[] cell_neighbors; delete[] num_faces_corners;
			delete[] faces_corners; delete[] faces_area;

			cell_neighbors = new_cell_neighbors;
			num_faces_corners = new_num_faces_corners;
			faces_corners = new_faces_corners;
			faces_area = new_faces_area;
			cell_neighbors_cap = new_cap;
			#pragma endregion
		}
	}

	delete[] additional_neighbors; delete[] additional_neighbors_x;
	delete[] spoke_dir; delete[] xface; delete[] seeds;

	//save_Voronoi_cell_ply("cell.ply", num_cell_neighbors, num_faces_corners, faces_corners);

	return 0;
	#pragma endregion
}

int MeshingSmartTree::get_3d_Voronoi_face(size_t iseed, size_t jseed, double* xface,
	                               size_t &num_face_corners, double**& face_corners,
	                               double &face_area, double &face_height)
{
	#pragma region Construct 3d Voronoi Face:

	if (_num_dim != 3) return 1;

	num_face_corners = 0;

	size_t* seeds = new size_t[4];
	seeds[0] = iseed; seeds[1] = jseed;

	double* si = _points[iseed]; double* sj = _points[jseed];

	double** e_basis = new double* [3];
	for (size_t i = 0; i < 3; i++) e_basis[i] = new double[3];

	double* xedge = new double[3];
	double* spoke_dir = new double[3];

	// first basis from si to sj
	double norm(0.0);
	for (size_t idim = 0; idim < 3; idim++)
	{
		e_basis[0][idim] = sj[idim] - si[idim];
		norm += e_basis[0][idim] * e_basis[0][idim];
	}
	norm = sqrt(norm);
	for (size_t idim = 0; idim < 3; idim++) e_basis[0][idim] /= norm;

	if (true)
	{
		#pragma region Inspect for degenercy:
		double hseed(0.0);
		for (size_t idim = 0; idim < 3; idim++)
		{
			double dx = xface[idim] - si[idim];
			hseed += dx * dx;
		}
		hseed = sqrt(hseed);
		size_t num_sphere_neighbors(0), sphere_neighbors_cap(10);
		size_t* sphere_neighbors = new size_t[sphere_neighbors_cap];
		get_tree_points_in_sphere(xface, hseed, num_sphere_neighbors, sphere_neighbors, sphere_neighbors_cap);

		bool found(false); size_t kseed(_num_points);
		for (size_t i = 0; i < num_sphere_neighbors; i++)
		{
			if (sphere_neighbors[i] == iseed) continue;
			else if (sphere_neighbors[i] == jseed) found = true;
		}


		if (!found)
		{
			num_face_corners = 0; face_area = 0.0;
			vcm_cout << "Warning Invalid Voronoi Face Point, iseed = "<< iseed << " , jseed = " << jseed << "!!!" << std::endl;

			for (size_t i = 0; i < 3; i++) delete[] e_basis[i];
			delete[] e_basis; delete[] seeds; delete[] xedge; delete[] spoke_dir;
			delete[] sphere_neighbors;

			return 1;
		}

		if (num_sphere_neighbors > 2)
		{
			#pragma region xface is a point on a Voronoi Edge:
			bool degenerate_face(true);
			for (size_t i = 0; i < num_sphere_neighbors; i++)
			{
				if (sphere_neighbors[i] == iseed || sphere_neighbors[i] == jseed) continue;

				kseed = sphere_neighbors[i];

				double* vec = new double[3]; double* sk = _points[kseed];
				for (size_t idim = 0; idim < 3; idim++) vec[idim] = sk[idim] - si[idim];

				norm = 0.0;
				for (size_t idim = 0; idim < 3; idim++) norm += vec[idim] * vec[idim];
				norm = sqrt(norm);
				for (size_t idim = 0; idim < 3; idim++) vec[idim] /= norm;

				double dot(0.0);
				for (size_t idim = 0; idim < 3; idim++) dot += vec[idim] * e_basis[0][idim];

				norm = 0.0;
				for (size_t idim = 0; idim < 3; idim++)
				{
					vec[idim] -= dot * e_basis[0][idim];;
					norm += vec[idim] * vec[idim];
				}
				norm = sqrt(norm);
				for (size_t idim = 0; idim < 3; idim++) e_basis[1][idim] = vec[idim] / norm;

				// third basis via cross product
				e_basis[2][0] = e_basis[0][1] * e_basis[1][2] - e_basis[0][2] * e_basis[1][1];
				e_basis[2][1] = e_basis[0][2] * e_basis[1][0] - e_basis[0][0] * e_basis[1][2];
				e_basis[2][2] = e_basis[0][0] * e_basis[1][1] - e_basis[0][1] * e_basis[1][0];

				delete[] vec;

				// throw two vectors along e1 and trim them
				bool degenerate_edge(true);
				for (size_t ii = 0; ii < 4; ii++)
				{
					if (ii >= 2 && degenerate_edge) break;

					if (ii == 0)
					{
						for (size_t idim = 0; idim < 3; idim++) spoke_dir[idim] = e_basis[2][idim];
					}
					else if (ii == 1)
					{
						for (size_t idim = 0; idim < 3; idim++) spoke_dir[idim] = - e_basis[2][idim];
					}
					else if (ii == 2)
					{
						for (size_t idim = 0; idim < 3; idim++) spoke_dir[idim] = e_basis[1][idim];
					}
					else
					{
						for (size_t idim = 0; idim < 3; idim++) spoke_dir[idim] = - e_basis[1][idim];
					}

					size_t neighbor_seed(_num_points); double spoke_length(DBL_MAX);

					if (ii < 2)
						trim_spoke(_points[iseed], 2, seeds, xface, spoke_dir, spoke_length, neighbor_seed);
					else
						trim_spoke(_points[iseed], 2, seeds, xedge, spoke_dir, spoke_length, neighbor_seed);

					if (spoke_length == DBL_MAX)
					{
						vcm_cout << "Unbounded Face!" << std::endl;

						face_area = DBL_MAX; face_height = 0.0;
						for (size_t idim = 0; idim < 3; idim++) face_height += (xface[idim] - si[idim]) * e_basis[0][idim];

						for (size_t i = 0; i < 3; i++) delete[] e_basis[i];
						delete[] e_basis; delete[] seeds; delete[] xedge; delete[] spoke_dir;
						delete[] sphere_neighbors;

						return 1;

					}

					if (spoke_length > 1E-10)
					{
						if (ii < 2)
						{
							for (size_t idim = 0; idim < 3; idim++) xedge[idim] = xface[idim] + 0.5 * spoke_length * spoke_dir[idim];
						}
						else
						{
							for (size_t idim = 0; idim < 3; idim++) xface[idim] = xedge[idim] + 0.5 * spoke_length * spoke_dir[idim];
						}


						if (ii < 2) degenerate_edge = false;
						else
						{
							degenerate_face = false;
							break;
						}
					}
				}
				if (!degenerate_face) break;
			}

			delete[] sphere_neighbors;

			if (degenerate_face)
			{
				num_face_corners = 0;
				face_area = 0.0; face_height = 0.0;
				for (size_t idim = 0; idim < 3; idim++) face_height += (xface[idim] - si[idim]) * e_basis[0][idim];

				for (size_t i = 0; i < 3; i++) delete[] e_basis[i];
				delete[] e_basis; delete[] seeds; delete[] xedge; delete[] spoke_dir;
				return 1;
			}
			#pragma endregion
		}
		else
		{
			#pragma region xface is a point on a Voronoi Face:
			delete[] sphere_neighbors;

			// Construct second basis via projection
			double best_norm(0.0);
			double* vec = new double[3];
			for (size_t i = 0; i < 3; i++)
			{
				for (size_t j = 0; j < 3; j++) vec[j] = 0.0;

				if (e_basis[0][i] > 0.0) vec[i] = 1.0;
				else                     vec[i] = -1.0;

				double dot(0.0);
				for (size_t idim = 0; idim < 3; idim++) dot += vec[idim] * e_basis[0][idim];

				norm = 0.0;
				for (size_t idim = 0; idim < 3; idim++)
				{
					vec[idim] -= dot * e_basis[0][idim];;
					norm += vec[idim] * vec[idim];
				}

				if (norm > best_norm)
				{
					best_norm = norm; norm = sqrt(norm);
					for (size_t idim = 0; idim < 3; idim++) e_basis[1][idim] = vec[idim] / norm;
				}
			}
			delete[] vec;

			// third basis via cross product
			e_basis[2][0] = e_basis[0][1] * e_basis[1][2] - e_basis[0][2] * e_basis[1][1];
			e_basis[2][1] = e_basis[0][2] * e_basis[1][0] - e_basis[0][0] * e_basis[1][2];
			e_basis[2][2] = e_basis[0][0] * e_basis[1][1] - e_basis[0][1] * e_basis[1][0];
			#pragma endregion
		}
		#pragma endregion
	}

	bool unbounded_face = false; bool duplicate_edges = false;

	size_t num_spokes(5);

	size_t num_joint_neighbors = 0;
	size_t joint_neighbors_cap(100);
	size_t* joint_neighbors = new size_t[joint_neighbors_cap];
	double* dart = new double[2];

	size_t jspoke(0);
	while (true)
	{
		size_t joint_neighbor_seed;
		if (jspoke < num_spokes || num_joint_neighbors == 1)
		{
			#pragma region get a joint neighbor through VoroSpokes:
			if (jspoke < num_spokes)
			{
				_rsampler.sample_uniformly_from_unit_sphere(dart, 2);
				for (size_t idim = 0; idim < 3; idim++)
				{
					spoke_dir[idim] = dart[0] * e_basis[1][idim] + dart[1] * e_basis[2][idim];
				}
			}
			else if (num_joint_neighbors == 1)
			{
				// flip last spoke
				for (size_t idim = 0; idim < 3; idim++) spoke_dir[idim] = -spoke_dir[idim];
			}

			jspoke++;

			double spoke_length(DBL_MAX);
			trim_spoke(_points[iseed], 2, seeds, xface, spoke_dir, spoke_length, joint_neighbor_seed);

			if (spoke_length == DBL_MAX)
			{
				// an unbounded face
				unbounded_face = true; break;
			}

			bool prior_joint_neighbor(false);
			for (size_t i = 0; i < num_joint_neighbors; i++)
			{
				if (joint_neighbors[i] == joint_neighbor_seed)
				{
					prior_joint_neighbor = true;
					break;
				}
			}

			if (prior_joint_neighbor)
			{
				if (jspoke <= num_spokes) continue; // We already collected that joint neighbor
				// A bug!
				vcm_cout << "Error::Invalid Spoke Trimming" << std::endl;
			}

			if (num_joint_neighbors == 1)
			{
				size_t neighbor_seed_index = joint_neighbors[0];
				double* q = _points[neighbor_seed_index];
				double dq_1(0.0), dq_2(0.0);
				for (size_t idim = 0; idim < 3; idim++)
				{
					double xmid = 0.5 * (si[idim] + sj[idim]);
					dq_1 += (q[idim] - xmid) * e_basis[1][idim];
					dq_2 += (q[idim] - xmid) * e_basis[2][idim];
				}
				double theta_o = _geom.get_point_angle(dq_1, dq_2);

				neighbor_seed_index = joint_neighbor_seed;
				q = _points[neighbor_seed_index];
				dq_1 = 0.0; dq_2 = 0.0;
				for (size_t idim = 0; idim < 3; idim++)
				{
					double xmid = 0.5 * (si[idim] + sj[idim]);
					dq_1 += (q[idim] - xmid) * e_basis[1][idim];
					dq_2 += (q[idim] - xmid) * e_basis[2][idim];
				}
				double theta = _geom.get_point_angle(dq_1, dq_2);

				if (fabs(theta - theta_o) < 1E-8) continue;
				if (fabs(theta - theta_o - 2 * PI) < 1E-8) continue;
			}

			joint_neighbors[num_joint_neighbors] = joint_neighbor_seed;
			num_joint_neighbors++;

			if (num_joint_neighbors == joint_neighbors_cap)
			{
				size_t* new_joint_neighbors = new size_t[joint_neighbors_cap * 2];
				for (size_t i = 0; i < num_joint_neighbors; i++) new_joint_neighbors[i] = joint_neighbors[i];
				delete[] joint_neighbors;
				joint_neighbors = new_joint_neighbors; joint_neighbors_cap *= 2;
			}
			#pragma endregion
		}
		else
		{
			#pragma region Construct Voronoi Face:

			double* theta = new double[num_joint_neighbors];
			double* theta_x = new double[num_joint_neighbors];
			double* theta_y = new double[num_joint_neighbors];

			for (size_t i = 0; i < num_joint_neighbors; i++)
			{
				#pragma region Retrieve Neibhors angles:
				size_t neighbor_seed_index = joint_neighbors[i];
				double* q = _points[neighbor_seed_index];
				double dq_1(0.0), dq_2(0.0);
				for (size_t idim = 0; idim < 3; idim++)
				{
					double xmid = 0.5 * (si[idim] + sj[idim]);
					dq_1 += (q[idim] - xmid) * e_basis[1][idim];
					dq_2 += (q[idim] - xmid) * e_basis[2][idim];
				}
				theta[i] = _geom.get_point_angle(dq_1, dq_2);
				theta_x[i] = dq_1; theta_y[i] = dq_2;
				#pragma endregion
			}

			// sort theta and theta_index
			quicksort(theta, theta_x, theta_y, joint_neighbors, 0, num_joint_neighbors - 1);

			// check that all angles are less than PI
			bool new_joint_neighbors_added(false);
			for (size_t i = 0; i < num_joint_neighbors; i++)
			{
				size_t im = num_joint_neighbors - 1; if (i > 0) im = i - 1;
				double alpha = theta[i] - theta[im];
				if (alpha < 0.0) alpha += 2 * PI;
				if (alpha > PI - 1E-5)
				{
					#pragma region Fix an invalid angle:
					double beta = theta[im] + 0.5 * alpha;
					if (beta > 2 * PI) beta -= 2 * PI;
					double vx = cos(beta); double vy = sin(beta);
					for (size_t idim = 0; idim < 3; idim++) spoke_dir[idim] = vx * e_basis[1][idim];
					for (size_t idim = 0; idim < 3; idim++) spoke_dir[idim] += vy * e_basis[2][idim];

					seeds[2] = joint_neighbors[im]; seeds[3] = joint_neighbors[i];

					double spoke_length(DBL_MAX);
					trim_spoke(_points[iseed], 4, seeds, xface, spoke_dir, spoke_length, joint_neighbor_seed);

					if (spoke_length == DBL_MAX)
					{
						// an infinit direction
						unbounded_face = true;
						break;
					}
					for (size_t idim = 0; idim < 3; idim++) xedge[idim] = xface[idim] + spoke_length * spoke_dir[idim];

					double dq_1(0.0), dq_2(0.0);
					for (size_t idim = 0; idim < 3; idim++)
					{
						double xmid = 0.5 * (si[idim] + sj[idim]);
						dq_1 += (_points[joint_neighbor_seed][idim] - xmid) * e_basis[1][idim];
						dq_2 += (_points[joint_neighbor_seed][idim] - xmid) * e_basis[2][idim];
					}

					double beta_seed = _geom.get_point_angle(dq_1, dq_2);
					double alpha_seed = beta_seed - theta[im];
					if (alpha_seed < 0) alpha_seed += 2 * PI;
					if (alpha_seed < 1E-10 || alpha_seed > alpha - 1E-10)
					{
						vcm_cout << "Error:: An invalid Voronoi Face Angle" << std::endl;
						break; // degenerate face can't collect neighbors via spokes
					}

					joint_neighbors[num_joint_neighbors] = joint_neighbor_seed;
					num_joint_neighbors++;
					new_joint_neighbors_added = true;
					if (num_joint_neighbors == joint_neighbors_cap)
					{
						size_t* new_joint_neighbors = new size_t[joint_neighbors_cap * 2];
						for (size_t i = 0; i < num_joint_neighbors; i++) new_joint_neighbors[i] = joint_neighbors[i];
						delete[] joint_neighbors;
						joint_neighbors = new_joint_neighbors; joint_neighbors_cap *= 2;
					}
					break;
					#pragma endregion
				}
			}

			if (unbounded_face)
			{
				delete[] theta; delete[] theta_x; delete[] theta_y;
				break;
			}

			if (new_joint_neighbors_added)
			{
				delete[] theta; delete[] theta_x; delete[] theta_y;
				continue;
			}

			// filter our redundant neighbors
			size_t ii(0);
			while (ii < num_joint_neighbors)
			{
				size_t im = num_joint_neighbors - 1; if (ii > 0) im = ii - 1;

				if (fabs(theta[ii] - theta[im]) < 1E-5 || fabs(fabs(theta[ii] - theta[im]) - 2 * PI) < 1E-5)
				{
					for (size_t i = ii + 1; i < num_joint_neighbors; i++)
					{
						theta[i - 1] = theta[i];
						theta_x[i - 1] = theta_x[i];
						theta_y[i - 1] = theta_y[i];
						joint_neighbors[i - 1] = joint_neighbors[i];
					}
					num_joint_neighbors--;
				}
				else ii++;
			}

			// We have a closed ring of neighbors
			num_face_corners = num_joint_neighbors;
			face_corners = new double*[num_face_corners];
			for (size_t i = 0; i < num_joint_neighbors; i++)
			{
				#pragma region Construct Face Corners and validate them:
				size_t im = num_joint_neighbors - 1; if (i > 0) im = i - 1;

				size_t im_index = joint_neighbors[im];
				size_t ip_index = joint_neighbors[i];

				face_corners[i] = new double[3];

				seeds[2] = im_index; seeds[3] = ip_index;
				get_Voronoi_vertex(4, seeds, face_corners[i]);

				double old_dst(0.0);
				for (size_t idim = 0; idim < 3; idim++)
				{
					spoke_dir[idim] = face_corners[i][idim] - xface[idim];
					old_dst += spoke_dir[idim] * spoke_dir[idim];
				}
				old_dst = sqrt(old_dst);

				for (size_t idim = 0; idim < 3; idim++) spoke_dir[idim] /= old_dst;

				// verify that this is indeed a Voronoi Vertex
				double spoke_length(DBL_MAX);
				trim_spoke(_points[iseed], 4, seeds, xface, spoke_dir, spoke_length, joint_neighbor_seed);

				if (true)
				{
					double* q = _points[joint_neighbor_seed];
					double dq_1(0.0), dq_2(0.0);
					for (size_t idim = 0; idim < 3; idim++)
					{
						double xmid = 0.5 * (si[idim] + sj[idim]);
						dq_1 += (q[idim] - xmid) * e_basis[1][idim];
						dq_2 += (q[idim] - xmid) * e_basis[2][idim];
					}
					double joint_neighbor_theta = _geom.get_point_angle(dq_1, dq_2);

					bool redundant(false);
					for (size_t ii = 0; ii < num_joint_neighbors; ii++)
					{
						if (fabs(theta[ii] - joint_neighbor_theta) < 1E-5 || fabs(fabs(theta[ii] - joint_neighbor_theta) - 2 * PI) < 1E-5)
						{
							redundant = true; break;
						}
					}
					if (redundant) continue;
				}

				if (spoke_length > old_dst - 1E-8) continue; // no trimming happend

				bool found(false);
				for (size_t ii = 0; ii < num_joint_neighbors; ii++)
				{
					if (joint_neighbors[ii] == joint_neighbor_seed)
					{
						found = true; break;
					}
				}

				if (!found)
				{
					#pragma region a joint neighbor is added
					joint_neighbors[num_joint_neighbors] = joint_neighbor_seed;
					num_joint_neighbors++;
					new_joint_neighbors_added = true;
					if (num_joint_neighbors == joint_neighbors_cap)
					{
						size_t* new_joint_neighbors = new size_t[joint_neighbors_cap * 2];
						for (size_t i = 0; i < num_joint_neighbors; i++) new_joint_neighbors[i] = joint_neighbors[i];
						delete[] joint_neighbors;
						joint_neighbors = new_joint_neighbors; joint_neighbors_cap *= 2;
					}
					for (size_t j = 0; j <= i; j++) delete[] face_corners[j];
					delete[] face_corners;
					break;
					#pragma endregion
				}
				else
				{
					// adjust location of corners without adding a new joint neighbor
					vcm_cout << "Warning in Voronoi Face Construction, iseed = " << iseed << ", jseed = " << jseed << "!!!!!" << std::endl;
					for (size_t j = 0; j <= i; j++) delete[] face_corners[j];
					delete[] face_corners;
					delete[] theta; delete[] theta_x; delete[] theta_y;
					for (size_t i = 0; i < 3; i++) delete[] e_basis[i];
					delete[] e_basis; delete[] seeds; delete[] xedge;
					delete[] dart; delete[] spoke_dir; delete[] joint_neighbors;

					num_face_corners = 0;
					face_area = 0.0;

					return 1;
					return get_3d_Voronoi_face(iseed, jseed, xface, num_face_corners, face_corners, face_area, face_height);
				}
				#pragma endregion
			}

			delete[] theta; delete[] theta_x; delete[] theta_y;
			if (new_joint_neighbors_added) continue;

			break;

			#pragma endregion
		}
	}

	if (unbounded_face)
	{
		face_area = DBL_MAX;
	}
	else if (num_face_corners == 0)
	{
		// degenerate face
		face_area = 0.0;
		for (size_t idim = 0; idim < 3; idim++) face_height += (xface[idim] - si[idim]) * e_basis[0][idim];
	}
	else
	{
		// save_Voronoi_face_ply("face.ply", num_joint_neighbors, face_corners);

		face_area = 0.0;
		for (size_t i = 1; i < num_face_corners - 1; i++)
		{
			#pragma region Estimating Face Area:
			size_t ip = i + 1;

			double ax = face_corners[i][0] - face_corners[0][0];
			double ay = face_corners[i][1] - face_corners[0][1];
			double az = face_corners[i][2] - face_corners[0][2];

			double bx = face_corners[ip][0] - face_corners[0][0];
			double by = face_corners[ip][1] - face_corners[0][1];
			double bz = face_corners[ip][2] - face_corners[0][2];

			double nx = ay * bz - az * by;
			double ny = az * bx - ax * bz;
			double nz = ax * by - ay * bx;

			double norm(sqrt(nx * nx + ny * ny + nz * nz));
			face_area += 0.5 * norm;
			#pragma endregion
		}

		face_height = 0.0;
		for (size_t idim = 0; idim < 3; idim++) face_height += (face_corners[0][idim] - si[idim]) * e_basis[0][idim];

		if (fabs(face_area) < 1E-10)
		{
			face_area = 0.0;
		}
	}

	for (size_t i = 0; i < 3; i++) delete[] e_basis[i];
	delete[] e_basis; delete[] seeds; delete[] xedge;
	delete[] dart; delete[] spoke_dir; delete[] joint_neighbors;

	return 0;
	#pragma endregion
}

bool MeshingSmartTree::get_Voronoi_vertex(size_t num_seeds, size_t* seeds, double* vertex)
{
	#pragma region Construct a Voronoi Vertex:
	if (num_seeds != _num_dim + 1) return false;

	double** A = new double*[_num_dim];
	double* b = new double[_num_dim];
	size_t iseed = seeds[0];
	for (size_t i = 1; i < num_seeds; i++)
	{
		size_t jseed = seeds[i];
		A[i - 1] = new double[_num_dim];
		for (size_t idim = 0; idim < _num_dim; idim++) A[i - 1][idim] = _points[jseed][idim] - _points[iseed][idim];

		b[i - 1] = 0.0;
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double xmid = 0.5 * (_points[jseed][idim] + _points[iseed][idim]);
			b[i - 1] += A[i - 1][idim] * xmid;
		}
	}
	_lin_solver.LS_QR_Solver(_num_dim, _num_dim, A, b, vertex);

	for (size_t i = 1; i < num_seeds; i++) delete[] A[i - 1];
	delete[] A; delete[] b;

	return true;
	#pragma endregion
}

int MeshingSmartTree::get_normal_component(size_t num_dim, size_t num_basis, double** basis, double* vect, double &norm)
{
	#pragma region Get Normal component to some basis:

	double* comp = new double[num_basis];

	// project point to current basis
	for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
	{
		comp[ibasis] = 0.0;
		for (size_t idim = 0; idim < num_dim; idim++) comp[ibasis] += vect[idim] * basis[ibasis][idim];
	}

	// get vector component orthogonal to current basis
	for (size_t ibasis = 0; ibasis < num_basis; ibasis++)
	{
		for (size_t idim = 0; idim < num_dim; idim++) vect[idim] -= comp[ibasis] * basis[ibasis][idim];
	}

	delete[] comp;

	norm = 0.0;
	for (size_t idim = 0; idim < num_dim; idim++) norm += vect[idim] * vect[idim];

	if (fabs(norm) < 1E-10) return 1;

	norm = 1.0 / sqrt(norm);
	for (size_t idim = 0; idim < num_dim; idim++) vect[idim] *= norm;

	return 0;
	#pragma endregion
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// kd-tree  Methods  /////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int MeshingSmartTree::kd_tree_balance_quicksort(size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted)
{
	#pragma region kd tree balance:
	kd_tree_quicksort_adjust_target_position(target_pos, left, right, active_dim, tree_nodes_sorted);

	// target position is correct .. add to tree
	if (_tree_origin == _capacity)    _tree_origin = tree_nodes_sorted[target_pos];
	else                              kd_tree_add_point(tree_nodes_sorted[target_pos]);

	/* recursion */
	active_dim++;
	if (active_dim == _num_dim) active_dim = 0;

	if (target_pos > left + 1)  kd_tree_balance_quicksort((left + target_pos - 1) / 2, left, target_pos - 1, active_dim, tree_nodes_sorted);
	else if (left < target_pos) kd_tree_add_point(tree_nodes_sorted[left]);

	if (target_pos + 1 < right)  kd_tree_balance_quicksort((target_pos + 1 + right) / 2, target_pos + 1, right, active_dim, tree_nodes_sorted);
	else if (right > target_pos) kd_tree_add_point(tree_nodes_sorted[right]);

	return 0;
	#pragma endregion
}

int MeshingSmartTree::kd_tree_quicksort_adjust_target_position(size_t target_pos, size_t left, size_t right, size_t active_dim, size_t* tree_nodes_sorted)
{
	#pragma region kd tree Quick sort pivot:
	size_t i = left, j = right;

	size_t pivot_seed = tree_nodes_sorted[(left + right) / 2];
	double pivot = _points[pivot_seed][active_dim];

	/* partition */
	while (i <= j)
	{
		while (_points[tree_nodes_sorted[i]][active_dim] < pivot)
			i++;
		while (_points[tree_nodes_sorted[j]][active_dim] > pivot)
			j--;
		if (i <= j)
		{
			size_t tmp_index = tree_nodes_sorted[i];
			tree_nodes_sorted[i] = tree_nodes_sorted[j];
			tree_nodes_sorted[j] = tmp_index;

			i++;
			if (j > 0) j--;
		}
	};

	/* recursion */
	if (j > 0 && left < j && left <= target_pos && j >= target_pos)
		kd_tree_quicksort_adjust_target_position(target_pos, left, j, active_dim, tree_nodes_sorted);
	if (i < right && i <= target_pos && right >= target_pos)
		kd_tree_quicksort_adjust_target_position(target_pos, i, right, active_dim, tree_nodes_sorted);

	return 0;
	#pragma endregion
}

int MeshingSmartTree::kd_tree_add_point(size_t seed_index)
{
	#pragma region kd tree add point:
	// insert sphere into tree
	size_t parent_index(_tree_origin); size_t d_index(0);
	size_t branch_height(1);
	while (true)
	{
		if (_points[seed_index][d_index] >  _points[parent_index][d_index])
		{
			if (_tree_right[parent_index] == parent_index)
			{
				_tree_right[parent_index] = seed_index;
				branch_height++;
				break;
			}
			else
			{
				parent_index = _tree_right[parent_index];
				branch_height++;
			}
		}
		else
		{
			if (_tree_left[parent_index] == parent_index)
			{
				_tree_left[parent_index] = seed_index;
				branch_height++;
				break;
			}
			else
			{
				parent_index = _tree_left[parent_index];
				branch_height++;
			}
		}
		d_index++;
		if (d_index == _num_dim) d_index = 0;
	}
	if (branch_height > _tree_max_height) _tree_max_height = branch_height;
	return 0;
	#pragma endregion
}

int MeshingSmartTree::kd_tree_get_seeds_in_sphere(double* x,                                                    // Sphere center
	                                       double r,                                                     // neighborhood radius
	                                       size_t d_index, size_t node_index,                            // indices to traverse the kd-tree
	                                       size_t &num_points_in_sphere, size_t* &points_in_sphere,      // number of points in sphere and their indices
	                                       size_t &capacity                                              // Size of points in sphere array
	                                      )
{
	#pragma region kd tree neighbor search:
	if (d_index == _num_dim) d_index = 0;

	if (!_marked_only || _marked[node_index])
	{
		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = _points[node_index][idim] - x[idim];
			dst_sq += dx * dx;
		}

		if (dst_sq < (1 + 2.0E-6) * r * r) add_entry(node_index, num_points_in_sphere, points_in_sphere, capacity);
	}

	bool check_right(false), check_left(false);
	double neighbor_min(x[d_index] - r), neighbor_max(x[d_index] + r);

	if (_tree_right[node_index] != node_index && neighbor_max >  _points[node_index][d_index]) check_right = true;
	if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;

	if (check_right) kd_tree_get_seeds_in_sphere(x, r, d_index + 1, _tree_right[node_index], num_points_in_sphere, points_in_sphere, capacity);
	if (check_left)  kd_tree_get_seeds_in_sphere(x, r, d_index + 1, _tree_left[node_index], num_points_in_sphere, points_in_sphere, capacity);

	return 0;
	#pragma endregion
}

int MeshingSmartTree::kd_tree_get_seeds_in_sphere(double* x,                                                    // Sphere center
	                                       double r,                                                     // neighborhood radius
	                                       size_t max_index,                                             // maximum index to be considered
	                                       size_t d_index, size_t node_index,                            // indices to traverse the kd-tree
	                                       size_t &num_points_in_sphere, size_t* &points_in_sphere,      // number of points in sphere and their indices
	                                       size_t &capacity                                              // Size of points in sphere array
)
{
	#pragma region kd tree neighbor search:
	if (d_index == _num_dim) d_index = 0;

	if (node_index <= max_index)
	{
		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = _points[node_index][idim] - x[idim];
			dst_sq += dx * dx;
		}
		if (dst_sq < r * r + DST_TOL) add_entry(node_index, num_points_in_sphere, points_in_sphere, capacity);
	}

	bool check_right(false), check_left(false);
	double neighbor_min(x[d_index] - r), neighbor_max(x[d_index] + r);

	if (_tree_right[node_index] != node_index && neighbor_max >  _points[node_index][d_index]) check_right = true;
	if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;

	if (check_right) kd_tree_get_seeds_in_sphere(x, r, max_index,  d_index + 1, _tree_right[node_index], num_points_in_sphere, points_in_sphere, capacity);
	if (check_left)  kd_tree_get_seeds_in_sphere(x, r, max_index, d_index + 1, _tree_left[node_index], num_points_in_sphere, points_in_sphere, capacity);

	return 0;
	#pragma endregion
}

int MeshingSmartTree::kd_tree_get_closest_seed(size_t seed_index,                                  // seed index
	                                    size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
	                                    size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
	                                   )
{
	#pragma region kd tree closest neighbor search:
	if (d_index == _num_dim) d_index = 0;

	if (!_marked_only || _marked[node_index])
	{
		if (seed_index != node_index)
		{
			double dst_sq(0.0);
			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				double dx = _points[seed_index][idim] - _points[node_index][idim];
				dst_sq += dx * dx;
			}

			if (dst_sq < closest_distance * closest_distance)
			{
				// add to neighbors
				closest_seed = node_index;
				closest_distance = sqrt(dst_sq);
			}
		}
	}

	bool check_right(false), check_left(false);

	double neighbor_min, neighbor_max;
	if (closest_distance == DBL_MAX)
	{
		if (_tree_right[node_index] != node_index) check_right = true;
		if (_tree_left[node_index] != node_index) check_left = true;
	}
	else
	{
		neighbor_min = _points[seed_index][d_index] - closest_distance; neighbor_max = _points[seed_index][d_index] + closest_distance;
		if (_tree_right[node_index] != node_index && neighbor_max >  _points[node_index][d_index]) check_right = true;
		if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;
	}


	if (check_right && check_left)
	{
		// check the half that x blongs to first
		if (_points[seed_index][d_index] > _points[node_index][d_index])
		{
			kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			neighbor_min = _points[seed_index][d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
			if (neighbor_min < _points[node_index][d_index])
			{
				kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
			}
		}
		else
		{
			kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
			neighbor_max = _points[seed_index][d_index] + closest_distance;
			if (neighbor_max > _points[node_index][d_index])
			{
				kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			}
		}
	}
	else if (check_right) kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
	else if (check_left)  kd_tree_get_closest_seed(seed_index, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);

	return 0;
	#pragma endregion
}

int MeshingSmartTree::kd_tree_get_closest_seed(double* x,                                          // point location
	                                    size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
	                                    size_t &closest_seed, double &closest_distance,      // index of closest seed and distance from it
	                                    size_t &num_nodes_visited
	                                   )
{
	#pragma region kd tree closest neighbor search:
	if (d_index == _num_dim) d_index = 0;

	if (!_marked_only || _marked[node_index])
	{
		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = x[idim] - _points[node_index][idim];
			dst_sq += dx * dx;
		}
		num_nodes_visited++;
		if (dst_sq < closest_distance * closest_distance)
		{
			// add to neighbors
			closest_seed = node_index;
			closest_distance = sqrt(dst_sq);
		}
	}

	bool check_right(false), check_left(false);

	double neighbor_min, neighbor_max;
	if (closest_distance == DBL_MAX)
	{
		if (_tree_right[node_index] != node_index) check_right = true;
		if (_tree_left[node_index] != node_index) check_left = true;
	}
	else
	{
		neighbor_min = x[d_index] - closest_distance; neighbor_max = x[d_index] + closest_distance;
		if (_tree_right[node_index] != node_index && neighbor_max > _points[node_index][d_index]) check_right = true;
		if (_tree_left[node_index] != node_index  && neighbor_min < _points[node_index][d_index]) check_left = true;
	}

	if (check_right && check_left)
	{
		// check the half that x blongs to first
		if (x[d_index] > _points[node_index][d_index])
		{
			kd_tree_get_closest_seed(x, d_index + 1, _tree_right[node_index], closest_seed, closest_distance, num_nodes_visited);
			neighbor_min = x[d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
			if (neighbor_min < _points[node_index][d_index])
			{
				kd_tree_get_closest_seed(x, d_index + 1, _tree_left[node_index], closest_seed, closest_distance, num_nodes_visited);
			}
		}
		else
		{
			kd_tree_get_closest_seed(x, d_index + 1, _tree_left[node_index], closest_seed, closest_distance, num_nodes_visited);
			neighbor_max = x[d_index] + closest_distance;
			if (neighbor_max > _points[node_index][d_index])
			{
				kd_tree_get_closest_seed(x, d_index + 1, _tree_right[node_index], closest_seed, closest_distance, num_nodes_visited);
			}
		}
	}
	else if (check_right) kd_tree_get_closest_seed(x, d_index + 1, _tree_right[node_index], closest_seed, closest_distance, num_nodes_visited);
	else if (check_left)  kd_tree_get_closest_seed(x, d_index + 1, _tree_left[node_index], closest_seed, closest_distance, num_nodes_visited);

	return 0;
	#pragma endregion
}

int MeshingSmartTree::kd_tree_get_closest_seed(double* x,                                          // point location
	                                    size_t max_index,                                   // maximum index to be considered
	                                    size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
	                                    size_t &closest_seed, double &closest_distance,      // index of closest seed and distance from it
	                                    size_t &num_nodes_visited
)
{
	#pragma region kd tree closest neighbor search:
	if (d_index == _num_dim) d_index = 0;

	if (node_index <= max_index)
	{
		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = x[idim] - _points[node_index][idim];
			dst_sq += dx * dx;
		}
		if (dst_sq < closest_distance * closest_distance)
		{
			// Update Closest Neighbor
			closest_seed = node_index;
			closest_distance = sqrt(dst_sq);
		}
	}

	num_nodes_visited++;

	bool check_right(false), check_left(false);

	double neighbor_min, neighbor_max;
	if (closest_distance == DBL_MAX)
	{
		if (_tree_right[node_index] != node_index) check_right = true;
		if (_tree_left[node_index] != node_index) check_left = true;
	}
	else
	{
		neighbor_min = x[d_index] - closest_distance; neighbor_max = x[d_index] + closest_distance;
		if (_tree_right[node_index] != node_index && neighbor_max > _points[node_index][d_index]) check_right = true;
		if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;
	}

	if (check_right && check_left)
	{
		// check the half that x blongs to first
		if (x[d_index] > _points[node_index][d_index])
		{
			kd_tree_get_closest_seed(x, max_index, d_index + 1, _tree_right[node_index], closest_seed, closest_distance, num_nodes_visited);
			neighbor_min = x[d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
			if (neighbor_min < _points[node_index][d_index])
			{
				kd_tree_get_closest_seed(x, max_index, d_index + 1, _tree_left[node_index], closest_seed, closest_distance, num_nodes_visited);
			}
		}
		else
		{
			kd_tree_get_closest_seed(x, max_index, d_index + 1, _tree_left[node_index], closest_seed, closest_distance, num_nodes_visited);
			neighbor_max = x[d_index] + closest_distance;
			if (neighbor_max > _points[node_index][d_index])
			{
				kd_tree_get_closest_seed(x, max_index, d_index + 1, _tree_right[node_index], closest_seed, closest_distance, num_nodes_visited);
			}
		}
	}
	else if (check_right) kd_tree_get_closest_seed(x, max_index, d_index + 1, _tree_right[node_index], closest_seed, closest_distance, num_nodes_visited);
	else if (check_left)  kd_tree_get_closest_seed(x, max_index, d_index + 1, _tree_left[node_index], closest_seed, closest_distance, num_nodes_visited);

	return 0;
	#pragma endregion
}

int MeshingSmartTree::kd_tree_get_closest_non_smooth_seed(bool directional_normal,                            // use directional point normal
	                                               double* x,                                          // point location
	                                               size_t num_normals, double** x_normals,             // point normals
	                                               double smoothness_threshold,                        // smoothness threshold
	                                               size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
	                                               size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
                                                  )
{
	#pragma region kd tree closest neighbor search:
	if (d_index == _num_dim) d_index = 0;

	bool smooth(false);
	double dot(0.0), dst(0.0);
	for (size_t idim = 0; idim < _num_dim; idim++)
	{
		double dx = x[idim] - _points[node_index][idim];
		dot += dx * _points_normal[node_index][idim];
		dst += dx * dx;
	}
	dst = sqrt(dst);
	if (dst > 1E-10)
	{
		dot /= dst;
		if (!directional_normal) dot = fabs(dot);
		else                     dot = sqrt(1.0 - dot * dot);
	}
	else
		dot = 1.0;

	if (dot > smoothness_threshold)
	{
		for (size_t i = 0; i < num_normals; i++)
		{
			// angle between two normals
			dot = 0.0;
			for (size_t idim = 0; idim < _num_dim; idim++) dot += x_normals[i][idim] * _points_normal[node_index][idim];
			if (!directional_normal)
				dot = fabs(dot);

			if (dot > smoothness_threshold)
			{
				smooth = true;
				break;
			}
		}
	}

	if (!smooth)
	{
		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = x[idim] - _points[node_index][idim];
			dst_sq += dx * dx;
		}

		if (dst_sq < closest_distance * closest_distance)
		{
			closest_seed = node_index;
			closest_distance = sqrt(dst_sq);
		}
	}

	bool check_right(false), check_left(false);

	double neighbor_min, neighbor_max;
	if (closest_distance == DBL_MAX)
	{
		if (_tree_right[node_index] != node_index) check_right = true;
		if (_tree_left[node_index] != node_index) check_left = true;
	}
	else
	{
		neighbor_min = x[d_index] - closest_distance; neighbor_max = x[d_index] + closest_distance;
		if (_tree_right[node_index] != node_index && neighbor_max > _points[node_index][d_index]) check_right = true;
		if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;
	}

	if (check_right && check_left)
	{
		// check the half that x blongs to first
		if (x[d_index] > _points[node_index][d_index])
		{
			kd_tree_get_closest_non_smooth_seed(directional_normal, x, num_normals, x_normals, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			neighbor_min = x[d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
			if (neighbor_min < _points[node_index][d_index])
			{
				kd_tree_get_closest_non_smooth_seed(directional_normal, x, num_normals, x_normals, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
			}
		}
		else
		{
			kd_tree_get_closest_non_smooth_seed(directional_normal, x, num_normals, x_normals, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);

			neighbor_max = x[d_index] + closest_distance;
			if (neighbor_max > _points[node_index][d_index])
			{
				kd_tree_get_closest_non_smooth_seed(directional_normal, x, num_normals, x_normals, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			}
		}
	}
	else if (check_right)
	{
		kd_tree_get_closest_non_smooth_seed(directional_normal, x, num_normals, x_normals, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
	}
	else if (check_left)
	{
		kd_tree_get_closest_non_smooth_seed(directional_normal, x, num_normals, x_normals, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
	}

	return 0;
	#pragma endregion
}

// used for estimating surface sizing function in VC_Wild
int MeshingSmartTree::kd_tree_get_closest_non_smooth_seed(double* x, double feature_TOL,                 // point location
	                                               double smoothness_threshold,                    // smoothness threshold
	                                               size_t d_index, size_t node_index,              // indices to traverse the kd-tree
	                                               size_t &closest_seed, double &closest_distance  // index of closest seed and distance from it
                                                  )
{
	#pragma region kd tree closest neighbor search:
	if (d_index == _num_dim) d_index = 0;

	bool smooth(false);
	double dot(0.0), dst(0.0);
	for (size_t idim = 0; idim < _num_dim; idim++)
	{
		double dx = x[idim] - _points[node_index][idim];
		dot += dx * _points_normal[node_index][idim];
		dst += dx * dx;
	}
	dst = sqrt(dst);
	if (dst > feature_TOL)
	{
		dot /= dst;
		dot = sqrt(1.0 - dot * dot);
	}
	else dot = 1.0;

	if (dot > smoothness_threshold) smooth = true;

	if (!smooth)
	{
		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = x[idim] - _points[node_index][idim];
			dst_sq += dx * dx;
		}

		if (dst_sq < closest_distance * closest_distance)
		{
			closest_seed = node_index;
			closest_distance = sqrt(dst_sq);
		}
	}

	bool check_right(false), check_left(false);

	double neighbor_min, neighbor_max;
	if (closest_distance == DBL_MAX)
	{
		if (_tree_right[node_index] != node_index) check_right = true;
		if (_tree_left[node_index] != node_index) check_left = true;
	}
	else
	{
		neighbor_min = x[d_index] - closest_distance; neighbor_max = x[d_index] + closest_distance;
		if (_tree_right[node_index] != node_index && neighbor_max > _points[node_index][d_index]) check_right = true;
		if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;
	}

	if (check_right && check_left)
	{
		// check the half that x blongs to first
		if (x[d_index] > _points[node_index][d_index])
		{
			kd_tree_get_closest_non_smooth_seed(x, feature_TOL, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			neighbor_min = x[d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
			if (neighbor_min < _points[node_index][d_index])
			{
				kd_tree_get_closest_non_smooth_seed(x, feature_TOL, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
			}
		}
		else
		{
			kd_tree_get_closest_non_smooth_seed(x, feature_TOL, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);

			neighbor_max = x[d_index] + closest_distance;
			if (neighbor_max > _points[node_index][d_index])
			{
				kd_tree_get_closest_non_smooth_seed(x, feature_TOL, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			}
		}
	}
	else if (check_right)
	{
		kd_tree_get_closest_non_smooth_seed(x, feature_TOL, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
	}
	else if (check_left)
	{
		kd_tree_get_closest_non_smooth_seed(x, feature_TOL, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
	}
	return 0;
	#pragma endregion
}


// used for estimating surface sizing function in VC_Wild
int MeshingSmartTree::kd_tree_get_closest_non_smooth_seed(double* x, double* x_normal,                 // point location
	                                               double smoothness_threshold,                    // smoothness threshold
	                                               size_t d_index, size_t node_index,              // indices to traverse the kd-tree
	                                               size_t& closest_seed, double& closest_distance  // index of closest seed and distance from it
                                                  )
{
	#pragma region kd tree closest neighbor search:
	if (d_index == _num_dim) d_index = 0;

	bool smooth(false);
	double dot_p(0.0), dot_x(0.0), dst(0.0);
	for (size_t idim = 0; idim < _num_dim; idim++)
	{
		double dx = x[idim] - _points[node_index][idim];
		dot_p += dx * _points_normal[node_index][idim];
		dot_x += dx * x_normal[idim];
		dst += dx * dx;
	}
	dst = sqrt(dst);
	if (dst > 1E-10)
	{
		dot_p /= dst;
		dot_p = sqrt(1.0 - dot_p * dot_p);
		dot_x /= dst;
		dot_x = sqrt(1.0 - dot_x * dot_x);
		if (dot_p > smoothness_threshold && dot_x > smoothness_threshold) smooth = true;
	}
	else smooth = true;

	if (!smooth)
	{
		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = x[idim] - _points[node_index][idim];
			dst_sq += dx * dx;
		}

		if (dst_sq < closest_distance * closest_distance)
		{
			closest_seed = node_index;
			closest_distance = sqrt(dst_sq);
		}
	}

	bool check_right(false), check_left(false);

	double neighbor_min, neighbor_max;
	if (closest_distance == DBL_MAX)
	{
		if (_tree_right[node_index] != node_index) check_right = true;
		if (_tree_left[node_index] != node_index) check_left = true;
	}
	else
	{
		neighbor_min = x[d_index] - closest_distance; neighbor_max = x[d_index] + closest_distance;
		if (_tree_right[node_index] != node_index && neighbor_max > _points[node_index][d_index]) check_right = true;
		if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;
	}

	if (check_right && check_left)
	{
		// check the half that x blongs to first
		if (x[d_index] > _points[node_index][d_index])
		{
			kd_tree_get_closest_non_smooth_seed(x, x_normal, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			neighbor_min = x[d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
			if (neighbor_min < _points[node_index][d_index])
			{
				kd_tree_get_closest_non_smooth_seed(x, x_normal, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
			}
		}
		else
		{
			kd_tree_get_closest_non_smooth_seed(x, x_normal, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);

			neighbor_max = x[d_index] + closest_distance;
			if (neighbor_max > _points[node_index][d_index])
			{
				kd_tree_get_closest_non_smooth_seed(x, x_normal, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			}
		}
	}
	else if (check_right)
	{
		kd_tree_get_closest_non_smooth_seed(x, x_normal, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
	}
	else if (check_left)
	{
		kd_tree_get_closest_non_smooth_seed(x, x_normal, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
	}
	return 0;
	#pragma endregion
}


// used for detecting sharp corners/edges in VC_Wild
int MeshingSmartTree::kd_tree_get_closest_non_smooth_seed(double* x, size_t num_normals, double** normals,                      // point location and normal
	                                               double smoothness_threshold,                                          // smoothness threshold
	                                               size_t d_index, size_t node_index,                                    // indices to traverse the kd-tree
	                                               size_t &closest_seed, double &closest_distance                        // index of closest seed and distance from it
                                                   )
{
	#pragma region kd tree closest neighbor search:
	if (d_index == _num_dim) d_index = 0;

	bool smooth(false);
	double dot(0.0), dst(0.0);
	for (size_t idim = 0; idim < _num_dim; idim++)
	{
		double dx = x[idim] - _points[node_index][idim];
		dot += dx * _points_normal[node_index][idim];
		dst += dx * dx;
	}
	if (dst > 1E-10)
	{
		dst = sqrt(dst); dot /= dst;
		dot = sqrt(1.0 - dot * dot);
	}
	else dot = 1.0;

	if (dot > smoothness_threshold)
	{
		for (size_t i = 0; i < num_normals; i++)
		{
			dot = 0.0;
			for (size_t idim = 0; idim < _num_dim; idim++) dot += normals[i][idim] * _points_normal[node_index][idim];
			dot = fabs(dot);
			if (dot > smoothness_threshold)
			{
				smooth = true;
				break;
			}
		}
	}

	if (!smooth)
	{
		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = x[idim] - _points[node_index][idim];
			dst_sq += dx * dx;
		}

		if (dst_sq < closest_distance * closest_distance)
		{
			closest_seed = node_index;
			closest_distance = sqrt(dst_sq);
		}
	}

	bool check_right(false), check_left(false);

	double neighbor_min, neighbor_max;
	if (closest_distance == DBL_MAX)
	{
		if (_tree_right[node_index] != node_index) check_right = true;
		if (_tree_left[node_index] != node_index) check_left = true;
	}
	else
	{
		neighbor_min = x[d_index] - closest_distance; neighbor_max = x[d_index] + closest_distance;
		if (_tree_right[node_index] != node_index && neighbor_max > _points[node_index][d_index]) check_right = true;
		if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;
	}

	if (check_right && check_left)
	{
		// check the half that x blongs to first
		if (x[d_index] > _points[node_index][d_index])
		{
			kd_tree_get_closest_non_smooth_seed(x, num_normals, normals, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			neighbor_min = x[d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
			if (neighbor_min < _points[node_index][d_index])
			{
				kd_tree_get_closest_non_smooth_seed(x, num_normals, normals, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
			}
		}
		else
		{
			kd_tree_get_closest_non_smooth_seed(x, num_normals, normals, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);

			neighbor_max = x[d_index] + closest_distance;
			if (neighbor_max > _points[node_index][d_index])
			{
				kd_tree_get_closest_non_smooth_seed(x, num_normals, normals, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			}
		}
	}
	else if (check_right)
	{
		kd_tree_get_closest_non_smooth_seed(x, num_normals, normals, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
	}
	else if (check_left)
	{
		kd_tree_get_closest_non_smooth_seed(x, num_normals, normals, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
	}
	return 0;
	#pragma endregion
}

// used for estimating curves sizing function in VC_Wild
int MeshingSmartTree::kd_tree_get_closest_non_smooth_edge_seed(double* x, double feature_TOL,                  // point location and feature tolerance
	                                                    double smoothness_threshold,                    // smoothness threshold
	                                                    size_t d_index, size_t node_index,              // indices to traverse the kd-tree
	                                                    size_t &closest_seed, double &closest_distance  // index of closest seed and distance from it
                                                       )
{
	#pragma region kd tree closest neighbor search:
	if (d_index == _num_dim) d_index = 0;

	bool smooth(false);
	double dot(0.0), dst(0.0);
	for (size_t idim = 0; idim < _num_dim; idim++)
	{
		double dx = x[idim] - _points[node_index][idim];
		dot += dx * _points_normal[node_index][idim];
		dst += dx * dx;
	}
	dst = sqrt(dst);

	if (dst > feature_TOL)
	{
		dot /= dst;
		dot = fabs(dot);
	}
	else dot = 1.0;

	if (dot > smoothness_threshold) smooth = true;

	if (!smooth)
	{
		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = x[idim] - _points[node_index][idim];
			dst_sq += dx * dx;
		}

		if (dst_sq < closest_distance * closest_distance)
		{
			closest_seed = node_index;
			closest_distance = sqrt(dst_sq);
		}
	}

	bool check_right(false), check_left(false);

	double neighbor_min, neighbor_max;
	if (closest_distance == DBL_MAX)
	{
		if (_tree_right[node_index] != node_index) check_right = true;
		if (_tree_left[node_index] != node_index) check_left = true;
	}
	else
	{
		neighbor_min = x[d_index] - closest_distance; neighbor_max = x[d_index] + closest_distance;
		if (_tree_right[node_index] != node_index && neighbor_max > _points[node_index][d_index]) check_right = true;
		if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;
	}

	if (check_right && check_left)
	{
		// check the half that x blongs to first
		if (x[d_index] > _points[node_index][d_index])
		{
			kd_tree_get_closest_non_smooth_edge_seed(x, feature_TOL, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			neighbor_min = x[d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
			if (neighbor_min < _points[node_index][d_index])
			{
				kd_tree_get_closest_non_smooth_edge_seed(x, feature_TOL, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
			}
		}
		else
		{
			kd_tree_get_closest_non_smooth_edge_seed(x, feature_TOL, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);

			neighbor_max = x[d_index] + closest_distance;
			if (neighbor_max > _points[node_index][d_index])
			{
				kd_tree_get_closest_non_smooth_edge_seed(x, feature_TOL, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			}
		}
	}
	else if (check_right)
	{
		kd_tree_get_closest_non_smooth_edge_seed(x, feature_TOL, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
	}
	else if (check_left)
	{
		kd_tree_get_closest_non_smooth_edge_seed(x, feature_TOL, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
	}
	return 0;
	#pragma endregion
}

// used for detecting sharp corners due to edges in VC_Wild
int MeshingSmartTree::kd_tree_get_closest_non_smooth_edge_seed(double* x, double* tangent,                      // point location and tangent
	                                                    double smoothness_threshold,                    // smoothness threshold
	                                                    size_t d_index, size_t node_index,              // indices to traverse the kd-tree
	                                                    size_t &closest_seed, double &closest_distance  // index of closest seed and distance from it
                                                       )
{
	#pragma region kd tree closest neighbor search:
	if (d_index == _num_dim) d_index = 0;

	bool smooth(false);
	double dot(0.0), dst(0.0);
	for (size_t idim = 0; idim < _num_dim; idim++)
	{
		double dx = x[idim] - _points[node_index][idim];
		dot += dx * _points_normal[node_index][idim];
		dst += dx * dx;
	}
	if (dst > 1E-10)
	{
		dst = sqrt(dst); dot /= dst;
		dot = fabs(dot);
	}
	else dot = 1.0;

	if (dot > smoothness_threshold)
	{
		dot = 0.0;
		for (size_t idim = 0; idim < _num_dim; idim++) dot += tangent[idim] * _points_normal[node_index][idim];
		dot = fabs(dot);
		if (dot > smoothness_threshold) smooth = true;
	}

	if (!smooth)
	{
		double dst_sq(0.0);
		for (size_t idim = 0; idim < _num_dim; idim++)
		{
			double dx = x[idim] - _points[node_index][idim];
			dst_sq += dx * dx;
		}

		if (dst_sq < closest_distance * closest_distance)
		{
			closest_seed = node_index;
			closest_distance = sqrt(dst_sq);
		}
	}

	bool check_right(false), check_left(false);

	double neighbor_min, neighbor_max;
	if (closest_distance == DBL_MAX)
	{
		if (_tree_right[node_index] != node_index) check_right = true;
		if (_tree_left[node_index] != node_index) check_left = true;
	}
	else
	{
		neighbor_min = x[d_index] - closest_distance; neighbor_max = x[d_index] + closest_distance;
		if (_tree_right[node_index] != node_index && neighbor_max > _points[node_index][d_index]) check_right = true;
		if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;
	}

	if (check_right && check_left)
	{
		// check the half that x blongs to first
		if (x[d_index] > _points[node_index][d_index])
		{
			kd_tree_get_closest_non_smooth_edge_seed(x, tangent, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			neighbor_min = x[d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
			if (neighbor_min < _points[node_index][d_index])
			{
				kd_tree_get_closest_non_smooth_edge_seed(x, tangent, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
			}
		}
		else
		{
			kd_tree_get_closest_non_smooth_edge_seed(x, tangent, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);

			neighbor_max = x[d_index] + closest_distance;
			if (neighbor_max > _points[node_index][d_index])
			{
				kd_tree_get_closest_non_smooth_edge_seed(x, tangent, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
			}
		}
	}
	else if (check_right)
	{
		kd_tree_get_closest_non_smooth_edge_seed(x, tangent, smoothness_threshold, d_index + 1, _tree_right[node_index], closest_seed, closest_distance);
	}
	else if (check_left)
	{
		kd_tree_get_closest_non_smooth_edge_seed(x, tangent, smoothness_threshold, d_index + 1, _tree_left[node_index], closest_seed, closest_distance);
	}
	return 0;
	#pragma endregion
}

int MeshingSmartTree::kd_tree_get_closest_seed(double* x,                                          // a point in space
	                                    size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
	                                    size_t num_exculded_seeds, size_t* exculded_seeds,  // Number of excluded seeds and their indices
	                                    size_t &closest_seed, double &closest_distance      // index of closest seed and distance from it
	                                   )
{
	#pragma region kd tree closest neighbor search:
	if (d_index == _num_dim) d_index = 0;

	if (!_marked_only || _marked[node_index])
	{
		if (!find_brute(node_index, exculded_seeds, num_exculded_seeds))
		{
			double dst_sq(0.0);
			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				double dx = x[idim] - _points[node_index][idim];
				dst_sq += dx * dx;
			}

			if (dst_sq < closest_distance * closest_distance)
			{
				// add to neighbors
				closest_seed = node_index;
				closest_distance = sqrt(dst_sq);
			}
		}
	}

	bool check_right(false), check_left(false);

	double neighbor_min, neighbor_max;
	if (closest_distance == DBL_MAX)
	{
		if (_tree_right[node_index] != node_index) check_right = true;
		if (_tree_left[node_index] != node_index) check_left = true;
	}
	else
	{
		neighbor_min = x[d_index] - closest_distance; neighbor_max = x[d_index] + closest_distance;
		if (_tree_right[node_index] != node_index && neighbor_max > _points[node_index][d_index]) check_right = true;
		if (_tree_left[node_index] != node_index && neighbor_min < _points[node_index][d_index]) check_left = true;
	}

	if (check_right && check_left)
	{
		// check the half that x blongs to first
		if (x[d_index] > _points[node_index][d_index])
		{
			kd_tree_get_closest_seed(x, d_index + 1, _tree_right[node_index], num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);
			neighbor_min = x[d_index] - closest_distance; // update neighbor min since closest_distance might have been shrunk
			if (neighbor_min < _points[node_index][d_index])
			{
				kd_tree_get_closest_seed(x, d_index + 1, _tree_left[node_index], num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);
			}

		}
		else
		{
			kd_tree_get_closest_seed(x, d_index + 1, _tree_left[node_index], num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);
			neighbor_max = x[d_index] + closest_distance;
			if (neighbor_max > _points[node_index][d_index])
			{
				kd_tree_get_closest_seed(x, d_index + 1, _tree_right[node_index], num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);
			}
		}
	}
	else if (check_right) kd_tree_get_closest_seed(x, d_index + 1, _tree_right[node_index], num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);
	else if (check_left)  kd_tree_get_closest_seed(x, d_index + 1, _tree_left[node_index], num_exculded_seeds, exculded_seeds, closest_seed, closest_distance);

	return 0;
	#pragma endregion
}


int MeshingSmartTree::kd_tree_trim_spoke(double* cell_seed, // could be a new seed(no in _points)
	                              size_t num_involved_seeds, size_t* involved_seeds, double* spoke_origin, double* spoke_dir, double origin_dst,
	                              size_t d_index, size_t node_index,                  // indices to traverse the kd-tree
	                              double &spoke_length, size_t &neighbor_seed         // spoke length and last trimming neighbor
                                 )
{
	#pragma region kd tree trim spoke search:
	if (d_index == _num_dim) d_index = 0;

	if (!_marked_only || _marked[node_index])
	{
		bool found(false);
		for (size_t i = 0; i < num_involved_seeds; i++)
		{
			if (node_index != involved_seeds[i]) continue;
			found = true; break;
		}

		if (!found)
		{
			double dot_o(0.0), dot_e(0.0);

			for (size_t idim = 0; idim < _num_dim; idim++)
			{
				double xmid = 0.5 * (_points[node_index][idim] + cell_seed[idim]);
				double normal = _points[node_index][idim] - cell_seed[idim];
				dot_o += (xmid - spoke_origin[idim]) * normal;
				dot_e += spoke_dir[idim] * normal;
			}

			if (dot_o < 0.0) dot_o = 0.0; // numerical errors since spoke origin belongs to the cell of spook seed

			if (fabs(dot_e) > 1E-8 && dot_e > 0.0)
			{
				double t = dot_o / dot_e;
				if (t > -1E-8 && t < spoke_length)
				{
					spoke_length = t; neighbor_seed = node_index;
				}
			}
		}
	}

	bool check_right(false), check_left(false);

	double H = fabs(_points[node_index][d_index] - spoke_origin[d_index]);

	if (H < 2 * (spoke_length + origin_dst) + 1E-8)
	{
		// search both branches
		if (_tree_right[node_index] != node_index) check_right = true;
		if (_tree_left[node_index] != node_index) check_left = true;
	}
	else
	{
		if (_tree_right[node_index] != node_index && spoke_origin[d_index] > _points[node_index][d_index]) check_right = true;
		if (_tree_left[node_index] != node_index && spoke_origin[d_index] < _points[node_index][d_index]) check_left = true;
	}

	if (check_right && check_left)
	{
		// check the half that x blongs to first
		if (spoke_origin[d_index] > _points[node_index][d_index])
		{
			kd_tree_trim_spoke(cell_seed, num_involved_seeds, involved_seeds, spoke_origin, spoke_dir, origin_dst, d_index + 1, _tree_right[node_index], spoke_length, neighbor_seed);
			if (H < 2 * (spoke_length + origin_dst) + 1E-10)
			{
				kd_tree_trim_spoke(cell_seed, num_involved_seeds, involved_seeds, spoke_origin, spoke_dir, origin_dst, d_index + 1, _tree_left[node_index], spoke_length, neighbor_seed);
			}
		}
		else
		{
			kd_tree_trim_spoke(cell_seed, num_involved_seeds, involved_seeds, spoke_origin, spoke_dir, origin_dst, d_index + 1, _tree_left[node_index], spoke_length, neighbor_seed);
			if (H < 2 * (spoke_length + origin_dst) + 1E-10)
			{
				kd_tree_trim_spoke(cell_seed, num_involved_seeds, involved_seeds, spoke_origin, spoke_dir, origin_dst, d_index + 1, _tree_right[node_index], spoke_length, neighbor_seed);
			}
		}
	}
	else if (check_right) kd_tree_trim_spoke(cell_seed, num_involved_seeds, involved_seeds, spoke_origin, spoke_dir, origin_dst, d_index + 1, _tree_right[node_index], spoke_length, neighbor_seed);
	else if (check_left)  kd_tree_trim_spoke(cell_seed, num_involved_seeds, involved_seeds, spoke_origin, spoke_dir, origin_dst, d_index + 1, _tree_left[node_index], spoke_length, neighbor_seed);
	return 0;
	#pragma endregion
}

int MeshingSmartTree::add_entry(size_t entry, size_t &num_entries, size_t* &I, size_t &capacity)
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

int MeshingSmartTree::add_entry(double* entry, size_t &num_entries, double** &I, size_t &capacity)
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

bool MeshingSmartTree::find_brute(size_t entry, size_t* I, size_t num_entries)
{
	#pragma region find using Brutal search:
	for (size_t i = 0; i < num_entries; i++) if (I[i] == entry) return true;
	return false;
	#pragma endregion
}


int MeshingSmartTree::quicksort(double* x, double* y1, double* y2, size_t* I, size_t left, size_t right)
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
			double tmp = x[i]; x[i] = x[j]; x[j] = tmp;

			tmp = y1[i]; y1[i] = y1[j]; y1[j] = tmp;

			tmp = y2[i]; y2[i] = y2[j]; y2[j] = tmp;

			size_t tmpi = I[i]; I[i] = I[j]; I[j] = tmpi;

			i++;
			if (j > 0) j--;
		}
	};

	/* recursion */

	if (j > 0 && left < j)
		quicksort(x, y1, y2, I, left, j);
	if (i < right)
		quicksort(x, y1, y2, I, i, right);

	return 0;
	#pragma endregion
}
