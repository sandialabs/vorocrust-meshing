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
// MeshingExodusIO.cpp                                            Last modified (10/28/2021) //
///////////////////////////////////////////////////////////////////////////////////////////////
#include "MeshingExodusIO.h"

#ifdef USE_EXODUS

MeshingExodusFileWriter::MeshingExodusFileWriter()
{
        activeFilePath = std::string();
        exoid = -1;

	blockTypeMap.insert({"EX_NODAL",EX_NODAL});
	blockTypeMap.insert({"EX_NODE_BLOCK",EX_NODE_BLOCK});
	blockTypeMap.insert({"EX_NODE_SET",EX_NODE_SET});
	blockTypeMap.insert({"EX_EDGE_BLOCK",EX_EDGE_BLOCK});
	blockTypeMap.insert({"EX_EDGE_SET",EX_EDGE_SET});
	blockTypeMap.insert({"EX_FACE_BLOCK",EX_FACE_BLOCK });
	blockTypeMap.insert({"EX_FACE_SET",EX_FACE_SET});
	blockTypeMap.insert({"EX_ELEM_BLOCK",EX_ELEM_BLOCK });
	blockTypeMap.insert({"EX_ELEM_SET",EX_ELEM_SET});
	blockTypeMap.insert({"EX_SIDE_SET",EX_SIDE_SET});
	blockTypeMap.insert({"EX_ELEM_MAP",EX_ELEM_MAP });
	blockTypeMap.insert({"EX_NODE_MAP",EX_NODE_MAP});
	blockTypeMap.insert({"EX_EDGE_MAP",EX_EDGE_MAP});
	blockTypeMap.insert({"EX_FACE_MAP",EX_FACE_MAP});
	blockTypeMap.insert({"EX_GLOBAL",EX_GLOBAL });
	blockTypeMap.insert({"EX_COORDINATE",EX_COORDINATE});
	blockTypeMap.insert({"EX_ASSEMBLY ",EX_ASSEMBLY });
	blockTypeMap.insert({"EX_BLOB",EX_BLOB});
	blockTypeMap.insert({"EX_INVALID",EX_INVALID});
}

MeshingExodusFileWriter::~MeshingExodusFileWriter()
{
        //If the file writer has an open file
        //when destructed, attempt to close it.
        if(hasOpenFile())
                closeFile();
}

void MeshingExodusFileWriter::openFile(std::string path)
{
        /*
        Opens an Exodus II file for writing using the Exodus library's
        C API.

        path: The filepath where we want to store the Exodus file.
        */

        if(fileExists(path)){
                std::cout << "[ExodusFileWriter] Warning, the file " <<
                        path << " already exists and will be overwritten."
                        << std::endl;
                //TODO: Deleting the old file for now. need to address an apparent file clobbering bug.
                std::remove(path.c_str());
        }

        const char *pathAsCString = path.c_str();
        int CPUWordSize = sizeof(double);
        int IOWordSize = sizeof(double);

        //The mode choices, explained...
        //EX_CLOBBER: Exodus will create a new file regardless of whether a file
        //            with that name already exists, overwriting the old file.
        //EX_64BIT_OFFSET: Allows Exodus to create a model that can store datasets
        //            larger than 2 gigabytes by putting it in 64-bit mode.
        int mode = EX_CLOBBER | EX_64BIT_OFFSET;

        int result = ex_create(pathAsCString,mode,&CPUWordSize,&IOWordSize);

        if(result < 0)
                throw std::runtime_error("[ExodusFileWriter] The Exodus library encountered an unexpected exception when creating the output file.");
        else {
                activeFilePath = path;
                exoid = result;
        }
}

void MeshingExodusFileWriter::initializeFile(std::string title, int numberOfDimensions,int numberOfNodes,
        int numberOfElements, int numberOfElementBlocks,
        int numberOfNodeSets, int numberOfSideSets,int numberOfFaces,int numberOfFaceBlocks){
        /* This function writes the initialization parameters to the exodus file.
        It must be called once (and only once!) before writing any data to the file.*/
        if(exoid < 0){
                throw std::runtime_error("[ExodusFileWriter] Called initializeFile without having an open file.");
        }
        /*int result = ex_put_init(exoid,title.c_str(),
                numberOfDimensions,numberOfNodes,numberOfElements,
                numberOfElementBlocks,numberOfNodeSets,numberOfSideSets);
        */
        ex_init_params par;
        strcpy(par.title, title.c_str());
        par.num_dim = numberOfDimensions;
        par.num_nodes = numberOfNodes;
        par.num_edge = 0;
        par.num_edge_blk = 0;
        par.num_face = numberOfFaces;
        par.num_face_blk = numberOfFaceBlocks;
        par.num_elem = numberOfElements;
        par.num_elem_blk = numberOfElementBlocks;
        par.num_node_sets = numberOfNodeSets;
        par.num_edge_sets = 0;
        par.num_face_sets = 0;
        par.num_side_sets = numberOfSideSets;
        par.num_elem_sets = 0;
        par.num_node_maps = 0;
        par.num_edge_maps = 0;
        par.num_face_maps = 0;
        par.num_elem_maps = 0;

        int result = ex_put_init_ext (exoid, &par);
        if(result < 0){
                throw std::runtime_error("[ExodusFileWriter] Called initializeFile but this routine has been called previously or the data file is read-only.");
        }
}

void MeshingExodusFileWriter::setCoordinatesOnAxis(std::list<double>* coordinates,int dimension){
        /* Sets values for nodal coordinates on the specified Axis.
        0 --> X
        1 --> Y
        2 --> Z

        Because the coordinates are floating point values, the application code must declare
        the arrays passed to be the appropriate type (float or double) to
        match the compute word size passed in ex_create() or ex_open().
        */
        if(coordinates == NULL)
                throw std::runtime_error("[ExodusFileWriter] Called setCoordinatesOnAxis, but the pointer passed must be non-NULL!");

        double coordinatesArray[coordinates->size()];
        std::copy(coordinates->begin(),coordinates->end(),coordinatesArray);

        int result;
        if(dimension == 0)
                result = ex_put_coord(exoid,coordinatesArray,NULL,NULL);
        else if(dimension == 1)
                result = ex_put_coord(exoid,NULL,coordinatesArray,NULL);
        else if(dimension == 2)
                result = ex_put_coord(exoid,NULL,NULL,coordinatesArray);
        else
                throw std::runtime_error("[ExodusFileWriter] Called setCoordinatesOnAxis, but specified an invalid dimension (0=X,1=Y,2=Z).");

        if(result < 0){
                throw std::runtime_error("[ExodusFileWriter] Called setCoordinates, but the underlying file is invalid, read-only, or not initalized properly.");
        }
}


void MeshingExodusFileWriter::setCoordinates(std::list<double>* xCoordinates, std::list<double>* yCoordinates, std::list<double>* zCoordinates){
        /* Sets values for nodal coordinates.

        Because the coordinates are floating point values, the application code must declare
        the arrays passed to be the appropriate type (float or double) to
        match the compute word size passed in ex_create() or ex_open().

        Parameters
        [in]	exoid	exodus file ID returned from a previous call to ex_create() or ex_open().
        [in]	x_coor	The X-coordinates of the nodes. If this is NULL, the X-coordinates will not be written.
        [in]	y_coor	The Y-coordinates of the nodes. These are stored only if num_dim > 1; otherwise, pass in NULL.
                        If this is NULL, the Y-coordinates will not be written.
        [in]	z_coor	The Z-coordinates of the nodes. These are stored only if num_dim > 2; otherwise, pass in NULL.
                        If this is NULL, the Z-coordinates will not be written
        */


        if(xCoordinates == NULL && yCoordinates == NULL && zCoordinates == NULL)
                throw std::runtime_error("[ExodusFileWriter] Called setCoordinates, but at least one pointer to a coordinate list must be non-NULL.");


        int xArray[(xCoordinates != NULL ? xCoordinates->size() : 0)];
        if (xCoordinates != NULL)
                std::copy(xCoordinates->begin(),xCoordinates->end(),xArray);

        int yArray[(yCoordinates != NULL ? yCoordinates->size() : 0)];
        if (yCoordinates != NULL)
                std::copy(yCoordinates->begin(),yCoordinates->end(),yArray);

        int zArray[(zCoordinates != NULL ? zCoordinates->size() : 0)];
        if (zCoordinates != NULL)
                std::copy(zCoordinates->begin(),zCoordinates->end(),zArray);

        int result = ex_put_coord(exoid,
                (xCoordinates != NULL ? xArray : NULL),
                (yCoordinates != NULL ? yArray : NULL),
                (zCoordinates != NULL ? zArray : NULL));
        if(result < 0){
                throw std::runtime_error("[ExodusFileWriter] Called setCoordinates, but the underlying file is invalid, read-only, or not initalized properly.");
        }

}

void MeshingExodusFileWriter::setCoordinateNames(std::list<std::string> coordinateNames){
        /* Sets the descriptions for the coordinate system.

        coordinateNames: A set of strings that
        describe the dimensions ike "X"/"Y"/"Z".
        */

        char** array = new char*[coordinateNames.size()];
        unsigned index = 0;
        for (std::list<std::string>::const_iterator it = coordinateNames.begin(); it != coordinateNames.end(); ++it) {
                array[index]= const_cast<char*>(it->c_str());
                index++;
        }

        int result = ex_put_coord_names(exoid, array);
        if(result < 0){
                throw std::runtime_error("[ExodusFileWriter] Called setCoordinateNames, but the underlying file is invalid, read-only, or not initalized properly.");
        }
}

void MeshingExodusFileWriter::defineBlock(std::string blockType, int blockIdentifier,std::string description,
        int numberOfEntries,int numberOfNodesPerEntry, int numberOfEdgesPerEntry,
        int numberOfFacesPerEntry, int numberOfAttributesPerEntry) {
        /* Defines an element/edge/face block in the Exodus file.*/

        ex_entity_type blockTypeEnum;
        auto search = blockTypeMap.find(blockType);
        if (search != blockTypeMap.end())
                blockTypeEnum = search->second;
        else
                throw std::runtime_error("[ExodusFileWriter] Called defineBlock with an unrecognized block type.");


        int result = ex_put_block(exoid,blockTypeEnum,blockIdentifier,description.c_str(),numberOfEntries,
        numberOfNodesPerEntry,numberOfEdgesPerEntry,numberOfFacesPerEntry,numberOfAttributesPerEntry);
        if(result < 0){
                throw std::runtime_error("[ExodusFileWriter] Called defineBlock, but the underlying file is invalid, read-only, or not initalized properly.");
        }

}

void MeshingExodusFileWriter::setBlockName(int blockIdentifier,std::string blockName) {

        int result = ex_put_name(exoid, EX_ELEM_BLOCK, blockIdentifier, const_cast<char*>(blockName.c_str()));
        if(result < 0){
                throw std::runtime_error("[ExodusFileWriter] Called setBlockName, but either the underlying file is invalid or the block id could not be found.");
        }
}


void MeshingExodusFileWriter::setBlockConnectivity(std::string blockType,int blockIdentifier,std::list<int>* nodeElementConnectivity,
        std::list<int>* elementEdgeConnectivity, std::list<int>* elementFaceConnectivity){
        /*
        Writes the connectivity array for a block.

        ex_put_conn() parameters

        exoid	exodus file id
        blk_type	type of block
        blk_id	id of block
        node_conn	node-element connectivity
        elem_edge_conn	element-edge connectivity (NULL if none)
        elem_face_conn	element-face connectivity (NULL if none)
        */

        ex_entity_type blockTypeEnum;
        auto search = blockTypeMap.find(blockType);
        if (search != blockTypeMap.end())
                blockTypeEnum = search->second;
        else
                throw std::runtime_error("[ExodusFileWriter] Called setBlockConnectivity with an unrecognized block type.");


        int numberOfEntriesInBlock;
        int numberOfNodesPerEntry;
        int numberOfEdgesPerEntry;
        int numberOfFacesPerEntry;
        int blockExists = ex_get_block(exoid,blockTypeEnum,blockIdentifier,NULL,&numberOfEntriesInBlock,&numberOfNodesPerEntry,&numberOfEdgesPerEntry,&numberOfFacesPerEntry,NULL);
        if(blockExists < 0)
                throw std::runtime_error("[ExodusFileWriter] Called setBlockConnectivity, no block exists with the given id.");
        int totalNumberOfNodes = numberOfEntriesInBlock * numberOfNodesPerEntry;
        if(nodeElementConnectivity != NULL && nodeElementConnectivity->size() > totalNumberOfNodes)
                throw std::runtime_error("[ExodusFileWriter] Called setBlockConnectivity, but there are more nodes listed in the node connectivity than there are nodes in the block!");
        int totalNumberOfEdges = numberOfEntriesInBlock * numberOfEdgesPerEntry;
        if(elementEdgeConnectivity != NULL && elementEdgeConnectivity->size() > totalNumberOfEdges)
                throw std::runtime_error("[ExodusFileWriter] Called setBlockConnectivity, but there are more edges listed in the edge connectivity than there are edges in the block!");
        int totalNumberOfFaces = numberOfEntriesInBlock * numberOfFacesPerEntry;
        if(elementFaceConnectivity != NULL && elementFaceConnectivity->size() > totalNumberOfFaces)
                throw std::runtime_error("[ExodusFileWriter] Called setBlockConnectivity, but there are more faces listed in the face connectivity than there are faces in the block!");


        int node_conn[(nodeElementConnectivity != NULL ? nodeElementConnectivity->size() : 0)];
        if (nodeElementConnectivity != NULL)
                std::copy(nodeElementConnectivity->begin(),nodeElementConnectivity->end(),node_conn);

        int elem_edge_conn[(elementEdgeConnectivity != NULL ? elementEdgeConnectivity->size() : 0)];
        if (elementEdgeConnectivity != NULL)
                std::copy(elementEdgeConnectivity->begin(),elementEdgeConnectivity->end(),elem_edge_conn);

        int elem_face_conn[(elementFaceConnectivity != NULL ? elementFaceConnectivity->size() : 0)];
        if (elementFaceConnectivity != NULL)
                std::copy(elementFaceConnectivity->begin(),elementFaceConnectivity->end(),elem_face_conn);

        int result = ex_put_conn(exoid,blockTypeEnum,blockIdentifier,
                (nodeElementConnectivity != NULL ? node_conn : NULL),
                (elementEdgeConnectivity != NULL ? elem_edge_conn : NULL),
                (elementFaceConnectivity != NULL ? elem_face_conn : NULL));

        if(result < 0){
                throw std::runtime_error("[ExodusFileWriter] Called setBlockConnectivity, but either the underlying file, block, or connectivity list is invalid.");
        }
}

void MeshingExodusFileWriter::setNodeMap(std::list<int>* nodeIdentifierMap){
         /* Writes out the node numbering map to the database; this allows
         the entity numbers to be non-contiguous.  This map is used for
         mapping between local and global node ids.

         nodeIdentifierMap: A list of integers indicating how
         local node ids should be mapped to global node ids.
         */ 

        if(!hasOpenFile())
                throw std::runtime_error("[ExodusFileWriter] Called setNodeMap, but the file writer has no file open for writing.");

        if(nodeIdentifierMap->size() == 0)
                throw std::runtime_error("[ExodusFileWriter] Called setNodeMap, but the list passed was empty.");

        int mapArray[nodeIdentifierMap->size()];
        std::copy(nodeIdentifierMap->begin(),nodeIdentifierMap->end(),mapArray);
        int result = ex_put_id_map(exoid,EX_NODE_MAP,&mapArray);
        if(result < 0){
                throw std::runtime_error("[ExodusFileWriter] Called setNodeMap, but the map is of a bad type, is empty, already exists, or the underlying library was otherwise unable to create or store the map.");
        }
}

void MeshingExodusFileWriter::setElementMap(std::list<int>* elementIdentifierMap){
          /* Writes out the element numbering map to the database; this allows
          the entity numbers to be non-contiguous.  This map is used for
          mapping between local and global element ids.

          elementIdentifierMap: A list of integers indicating how
          local element ids should be mapped to global node ids.
          */

          if(!hasOpenFile())
                  throw std::runtime_error("[ExodusFileWriter] Called setNodeMap, but the underlying file is not open for writing.");

          if(elementIdentifierMap->size() == 0)
                  throw std::runtime_error("[ExodusFileWriter] Called setElementMap, but the list passed was empty.");

          int mapArray[elementIdentifierMap->size()];
          std::copy(elementIdentifierMap->begin(),elementIdentifierMap->end(),mapArray);
          int result = ex_put_id_map(exoid,EX_ELEM_MAP,&mapArray);
          if(result < 0){
                  throw std::runtime_error("[ExodusFileWriter] Called setElementMap, but the map is of a bad type, is empty, already exists, or the underlying library was otherwise unable to create or store the map.");
          }
}


void MeshingExodusFileWriter::defineNodeSet(int nodeSetId,std::string nodeSetName, std::list<int>* nodesInSet){
        int numberOfNodesInSet = nodesInSet->size();
        if(numberOfNodesInSet == 0)
                throw std::runtime_error("[ExodusFileWriter] Called defineNodeSet, but the list of nodes is empty.");

        int nodeSetArray[nodesInSet->size()];
        std::copy(nodesInSet->begin(),nodesInSet->end(),nodeSetArray);

        ex_put_set_param(exoid, EX_NODE_SET, nodeSetId, nodesInSet->size(), 0);

        int createSetResult = ex_put_set(exoid, EX_NODE_SET, nodeSetId, nodeSetArray, NULL);
        if(createSetResult < 0)
                throw std::runtime_error("[ExodusFileWriter] Called defineNodeSet, but the underlying file is either invalid or improperly initialized.");

        int nameSetResult = ex_put_name(exoid, EX_NODE_SET, nodeSetId, const_cast<char*>(nodeSetName.c_str()));
        if(createSetResult < 0)
                throw std::runtime_error("[ExodusFileWriter] Called defineNodeSet, but the writer was unable to assign a name to the node set, which indicates that set creation failed.");
}


void MeshingExodusFileWriter::defineSideSet(int sideSetId,std::string sideSetName, std::list<int>* elementsInSet, std::list<int>* sidesInSet){
        int numberOfSidesInSet = sidesInSet->size();
        if(numberOfSidesInSet == 0)
                throw std::runtime_error("[ExodusFileWriter] Called defineSideSet, but the list of sides is empty.");

        int elementSetArray[elementsInSet->size()];
        std::copy(elementsInSet->begin(),elementsInSet->end(),elementSetArray);

        int sideSetArray[sidesInSet->size()];
        std::copy(sidesInSet->begin(),sidesInSet->end(),sideSetArray);

        ex_put_set_param(exoid, EX_SIDE_SET, sideSetId, elementsInSet->size(), 0);
        ex_put_set_param(exoid, EX_SIDE_SET, sideSetId, sidesInSet->size(), 0);


        int createSetResult = ex_put_set(exoid, EX_SIDE_SET, sideSetId, sideSetArray, elementSetArray);
        if(createSetResult < 0)
                throw std::runtime_error("[ExodusFileWriter] Called defineSideSet, but the underlying file is either invalid or improperly initialized.");

        int nameSetResult = ex_put_name(exoid, EX_SIDE_SET, sideSetId, const_cast<char*>(sideSetName.c_str()));
        if(createSetResult < 0)
                throw std::runtime_error("[ExodusFileWriter] Called defineSideSet, but the writer was unable to assign a name to the side set, which indicates that set creation failed.");
}

void MeshingExodusFileWriter::closeFile() {
        /* Checks whether the file writer has a file open via the Exodus library*/
        if(!hasOpenFile()){
                throw std::runtime_error("[ExodusFileWriter] Attempted to close out an Exodus file, but the writer doesn't have a file open!");
        }
        ex_close(exoid);
        exoid = -1;
        activeFilePath = std::string();
}

bool MeshingExodusFileWriter::hasOpenFile() {
        /* Checks whether the file writer has a file open via the Exodus library.*/
        return exoid >= 0;
}

std::string MeshingExodusFileWriter::getActiveFilePath(){
        return activeFilePath;
}

int MeshingExodusFileWriter::getExodusFileIdentifier() {
        /*
        Return the library-internal file identifier held by the ExodusFileWriter.
        */
        return exoid;
}


bool MeshingExodusFileWriter::fileExists(std::string path){
        /*
        Private convenience function that checks whether a file path already exists.

        path: The filepath that the caller wants to test.
        */
        std::ifstream infile(path);
        return infile.good();
}
#endif //USE_EXODUS
