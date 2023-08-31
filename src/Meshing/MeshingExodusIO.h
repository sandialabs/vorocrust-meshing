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
// MeshingExodusIO.h                                              Last modified (10/28/2021) //
///////////////////////////////////////////////////////////////////////////////////////////////
#include "Config.h"

#ifdef USE_EXODUS

#include "exodusII.h"
#include "MeshingCommon.h"

#include <cmath>
#include <vector>
#include <list>
#include <utility>
#include <stack>
#include <map>
#include <set>
#include <tuple>
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


class MeshingExodusFileWriter {

public:
        MeshingExodusFileWriter();
        ~MeshingExodusFileWriter();
        void openFile(std::string path);
        void initializeFile(std::string title, int numberOfDimensions,int numberOfNodes,
                int numberOfElements,int numberOfElementBlocks,
                int numberOfNodeSets, int numberOfSideSets,int numberOfFaces=0,int numberOfFaceBlocks=0);
        void setCoordinateNames(std::list<std::string> coordinateNames);
        void setCoordinates(std::list<double>* xCoordinates, std::list<double>* yCoordinates, std::list<double>* zCoordinates);
        void setCoordinatesOnAxis(std::list<double>* coordinates,int dimension);
        void defineBlock(std::string blockType, int blockIdentifier,std::string description,
        int numberOfEntries,int numberOfNodesPerEntry, int numberOfEdgesPerEntry,
        int numberOfFacesPerEntry, int numberOfAttributesPerEntry);
        void setBlockName(int blockIdentifier,std::string blockName);
        void setBlockConnectivity(std::string blockType,int blockIdentifier,std::list<int>* nodeElementConnectivity,
        std::list<int>* elementEdgeConnectivity, std::list<int>* elementFaceConnectivity);
        void setNodeMap(std::list<int>* nodeIdentifierMap);
        void setElementMap(std::list<int>* elementIdentifierMap);
        void defineNodeSet(int nodeSetId,std::string nodeSetName, std::list<int>* nodesInSet);
        void defineSideSet(int sideSetId,std::string sideSetName, std::list<int>* elementsInSet, std::list<int>* sidesInSet);
        void closeFile();
        bool hasOpenFile();
        std::string getActiveFilePath();
        int getExodusFileIdentifier();
        bool fileExists(std::string path);

private:
        std::string activeFilePath;
        int exoid;
        std::map<std::string, ex_entity_type> blockTypeMap;
};




#endif //USE_EXODUS
