#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
//  MeshingOptionParser.h                                          Last modified (04/15/2020)//
///////////////////////////////////////////////////////////////////////////////////////////////
#endif // DOXYGEN_SHOULD_SKIP_THIS

#ifndef _OPTIONPARSER_H_
#define _OPTIONPARSER_H_

#include <algorithm>
#include <exception>
#include <iostream>
#include <vector>


class MeshingOptionParser
{
    public:

	MeshingOptionParser(const int argc, char* argv[])
	: appname(argv[0]), argvec(argv+1, argv+argc)
	{
	}


	virtual ~MeshingOptionParser() = default;


	virtual bool parse_options(void)
	{
		if( option_exists("-arg") )
		{
			this->method   = 1;
			this->filename = get_option_value("-arg", true);
			return true;
		}
		return false;
	}


	virtual void PrettyPrint()
	{
		for(auto arg: this->argvec)
		{
			std::cout << "- " << arg << std::endl;
		}
	}


	virtual void usage() const {}


	void parse(void)
	{
		// Display help and exit if no arguments are provided or -h is given.
		if( 0==argvec.size() || option_exists("-h") )
		{
			usage();
			exit(0);
		}

		// Determine method & filename
		if( !parse_options() )
		{
			std::cout << "Error::Missing or invalid <method> parameter." << std::endl;
			usage();
			exit(1);
		}
	}


	const char* get_option_value(const std::string& option, bool throw_if_missing=false)
	{
		std::vector<std::string>::const_iterator it;
		it = std::find(argvec.begin(), argvec.end(), option);
		if(it != argvec.end() && ++it != argvec.end())
		{
			return it->c_str();
		}
		if(throw_if_missing)
		{
			throw std::runtime_error("ERROR: Missing required parameter to program argument " + option);
		}
		return NULL;
	}


	bool option_exists(const std::string& option)
	{
		return std::find(argvec.begin(), argvec.end(), option) != argvec.end();
	}


	//
	// Getters
	//
	const int get_method(void) const
	{
		return this->method;
	}


	const std::string get_filename(void) const
	{
		return this->filename;
	}

	//
	// Class Data
	//
	public:

	std::string appname;
	std::vector<std::string> argvec;

	int method;
	std::string filename;

};    // class MeshingOptionParser


#endif // _OPTIONPARSER_H_
