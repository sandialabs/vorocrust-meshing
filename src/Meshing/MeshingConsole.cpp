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
//  MeshingConsole.cpp                                             Last modified (07/11/2020)//
///////////////////////////////////////////////////////////////////////////////////////////////
#endif // DOXYGEN_SHOULD_SKIP_THIS

#include "Version.h"
#include "MeshingConsole.h"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <time.h>



/**
 * Returns a time_t object containing the current system time.
 * A time_t object contains a UNIX style datetime value that is consistent with
 * the number of seconds since epoch.
 *
 * @return time_t object containing the current time.
 */
time_t time_now_()
{
    return std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
}


/**
 * Returns the current time as a formatted string like:
 * - YYYY-MM-DD HH:MM:SS UTC_OFFSET (i.e., 2020-03-01 12:30:02 -600)
 * - YYYY-MM-DD_HH-MM-SSUTC_OFFSET (i.e., 2020-03-01_12-30-02-600)
 *
 * @param[in] time the UNIX time value to be converted.
 * @param[in] no_spaces_or_colons If true then format the string without whitespaces.
 * @return string containing the formatted time.
 */
std::string time_as_string_(time_t time, bool no_spaces_or_colons=false)
{
    std::stringstream ss;
    struct tm datetime;

    #ifdef _WIN32
    localtime_s(&datetime, &time);  // For Windows Systems
    #else
    localtime_r(&time, &datetime);  // For POSIX Systems
    #endif

    if(no_spaces_or_colons)
    {
        //ss << std::put_time(localtime(&time), "%F_%H-%M-%S%z");
        ss << std::put_time(&datetime, "%F_%H-%M-%S%z");
    }
    else
    {
        //ss << std::put_time(localtime(&time), "%F %H:%M:%S %z");
        ss << std::put_time(&datetime, "%F %H:%M:%S %z");
    }
    return ss.str();
}



MeshingConsole::MeshingConsole()
{
    this->fs.open("VoroCrustLog.txt", std::fstream::out);

    this->time_start = time_now_();

}


MeshingConsole::~MeshingConsole()
{
    this->fs.close();
}



std::string MeshingConsole::generate_banner() const
{
    std::stringstream ss;

    ss << "***********************************************************" << std::endl
       << "*                                                         *" << std::endl
       << "*  Program    : VoroCrust                                 *" << std::endl
       << "*  Version    : " << std::setw(40) << std::left << VOROCRUST_VERSION << "  *" << std::endl
       << "*  Git SHA1   : " << std::setw(40) << std::left << GIT_COMMIT_SHA1 << "  *" << std::endl;
    if("UNKNOWN" != GIT_COMMIT_DATE)
    {
       ss << "*  Commit Date: " << std::setw(40) << std::left << time_as_string_(std::atoi(GIT_COMMIT_DATE)) << "  *" << std::endl;
    }
    else
    {
        ss << "*  Commit Date: " << std::setw(40) << std::left << GIT_COMMIT_DATE << "  *" << std::endl;
    }
    ss << "*  Run Date   : " << std::setw(40) << std::left << time_as_string_(this->time_start) << "  *" << std::endl
       << "*  Website    : https://vorocrust.sandia.gov              *" << std::endl
       << "*                                                         *" << std::endl
       << "***********************************************************" << std::endl;

    return ss.str();
}

std::string MeshingConsole::get_current_datetime() const
{
    return time_as_string_( time_now_() );
}


MeshingConsole& MeshingConsole::operator()()
{
    std::stringstream ss;

    this->fs << ss.str();
    std::cout << ss.str();
    return *this;
}


MeshingConsole& MeshingConsole::operator<<(std::ostream& (*value)(std::ostream& o))
{
    this->fs << value;
    std::cout << value;
    return *this;
}

