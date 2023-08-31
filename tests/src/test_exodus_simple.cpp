#include "Config.h"
#include <stdexcept>

#ifdef USE_EXODUS

#include "MeshingExodusIO.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <exception>
#include <string>
#define NDEBUG

void test_ExodusFileWriter_isDirectlyConstructible(){
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
}

void test_ExodusFileWriter_isDirectlyDestructible(){
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        delete writer;
}

void test_ExodusFileWriter_DoesNotDetectNonExistentFile(){
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        assert(!writer->fileExists("nonexistent_file.exo"));
        delete writer;
}

void test_ExodusFileWriter_canOpenNewFile(){
        std::remove("temporary.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        assert(!writer->hasOpenFile());
        writer->openFile("temporary.exo");
        assert(writer->fileExists("temporary.exo"));
        assert(writer->hasOpenFile());
        assert(writer->getActiveFilePath().compare("temporary.exo") == 0);
        std::remove("temporary.exo");
        delete writer;
}

void test_ExodusFileWriter_canOpenAndCloseArbitrarilyManyFiles(){
        std::list<MeshingExodusFileWriter*> container;

        for(int i = 0; i < 50; i++) {
                MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
                std::string fileName = "filetest"+std::to_string(i)+".exo";
                writer->openFile(fileName);
                container.insert(container.end(), writer);
        }

        int i = 0;
        for (MeshingExodusFileWriter* writer : container) {
                writer->closeFile();
                delete writer;
                std::string fileName = "filetest"+std::to_string(i)+".exo";
                std::remove(const_cast<char*>(fileName.c_str()));
                i++;
        }

}


void test_ExodusFileWriter_canOpenExistingFile(){
        std::remove("existing_file.exo");
        std::ofstream f("existing_file.exo");
        f << "This data will be overwritten." << std::endl;
        f.close();
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        assert(writer->fileExists("existing_file.exo"));
        writer->openFile("existing_file.exo");
        assert(writer->hasOpenFile());
        assert(writer->getActiveFilePath().compare("existing_file.exo") == 0);
        std::remove("existing_file.exo");
}

void test_ExodusFileWriter_canCloseFile(){
        std::remove("file_to_be_closed.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("file_to_be_closed.exo");
        assert(writer->hasOpenFile());
        writer->closeFile();
        assert(!writer->hasOpenFile());
        assert(writer->getActiveFilePath().length() == 0);
        assert(writer->fileExists("file_to_be_closed.exo"));
        std::remove("file_to_be_closed.exo");
}

void test_ExodusFileWriter_closingNonexistentFileCausesException(){
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        try {
                writer->closeFile();
        } catch(std::exception& e){
                return;
        }
        assert(false && "The call in the try block should have generated an exception.");
}

void test_ExodusFileWriter_canInitializeFile(){
        std::remove("file_to_be_initialized.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("file_to_be_initialized.exo");
        int num_dim, num_nods, num_el, num_el_blk, num_ns, num_ss;

        num_dim = 3; num_nods = 46;
        num_el = 5; num_el_blk = 5;
        num_ns = 2; num_ss = 5;

        writer->initializeFile("Database Title", num_dim,
                num_nods, num_el, num_el_blk, num_ns, num_ss);

        writer->closeFile();
        std::remove("file_to_be_initialized.exo");
        delete writer;
}

void test_ExodusFileWriter_cantCallInitializeFileMoreThanOnce(){
        std::remove("twice_initialized.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("twice_initialized.exo");
        int num_dim, num_nods, num_el, num_el_blk, num_ns, num_ss;

        num_dim = 3; num_nods = 46;
        num_el = 5; num_el_blk = 5;
        num_ns = 2; num_ss = 5;

        writer->initializeFile("Database Title", num_dim,
        num_nods, num_el, num_el_blk, num_ns, num_ss);

        try {
                writer->initializeFile("Database Title", num_dim,
                num_nods, num_el, num_el_blk, num_ns, num_ss);
        } catch(std::exception& e){
                writer->closeFile();
                std::remove("twice_initialized.exo");
                delete writer;
                return;
        }
        writer->closeFile();
        std::remove("twice_initialized.exo");
        delete writer;
        assert(false && "The call in the try block should have generated an exception.");
}

void test_ExodusFileWriter_cantCallInitializeFileWithNoOpenFile(){
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        int num_dim, num_nods, num_el, num_el_blk, num_ns, num_ss;

        num_dim = 3; num_nods = 46;
        num_el = 5; num_el_blk = 5;
        num_ns = 2; num_ss = 5;

        try {
                writer->initializeFile("Database Title", num_dim,
                num_nods, num_el, num_el_blk, num_ns, num_ss);
        } catch(std::exception& e){
                return;
        }
        assert(false && "The call in the try block should have generated an exception.");
}

void test_ExodusFileWriter_canSetCoordinateNames(){
        std::remove("coordinateNameTest.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("coordinateNameTest.exo");
        int num_dim, num_nods, num_el, num_el_blk, num_ns, num_ss;

        num_dim = 3; num_nods = 46;
        num_el = 5; num_el_blk = 5;
        num_ns = 2; num_ss = 5;

        writer->initializeFile("Database Title", num_dim,
        num_nods, num_el, num_el_blk, num_ns, num_ss);


        std::list<std::string> coordinateNames = {"X", "Y", "Z"};

        writer->setCoordinateNames(coordinateNames);
        writer->closeFile();
        delete writer;
        std::remove("coordinateNameTest.exo");
}

void test_ExodusFileWriter_mustInitializeFileBeforeSettingCoordinateNames(){
        std::remove("coordinateNameFailure.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("coordinateNameFailure.exo");

        try {
                std::list<std::string> coordinateNames = {"X", "Y", "Z"};
                writer->setCoordinateNames(coordinateNames);
        } catch(std::exception& e){
                writer->closeFile();
                delete writer;
                std::remove("coordinateNameFailure.exo");
                return;
        }
        writer->closeFile();
        delete writer;
        std::remove("coordinateNameFailure.exo");
        assert(false && "The call in the try block should have generated an exception.");
}

void test_ExodusFileWriter_canSetCoordinates(){
        std::remove("coordinateTest.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("coordinateTest.exo");
        int numberOfDimensions = 3;
        int numberOfNodes = 4;
        int numberOfElements = 1;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 0;
        int numberOfSideSets = 0;

        writer->initializeFile("Coordinate Test", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks, numberOfNodeSets, numberOfSideSets);

        std::list<double> xCoordinates{0.0,1.0,0.0,1.0};
        std::list<double> yCoordinates{1.0,0.0,0.0,1.0};
        std::list<double> zCoordinates{0.0,0.0,0.0,0.0};

        writer->setCoordinatesOnAxis(&xCoordinates,0);
        writer->setCoordinatesOnAxis(&yCoordinates,1);
        writer->setCoordinatesOnAxis(&zCoordinates,2);

        writer->closeFile();
        std::remove("coordinateTest.exo");
}

void test_ExodusFileWriter_canCreateBlock(){
        std::remove("blockTest.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("blockTest.exo");
        int num_dim, num_nodes, num_el, num_el_blk, num_ns, num_ss;

        num_dim = 3; num_nodes = 33;
        num_el = 7; num_el_blk = 7;
        num_ns = 2; num_ss = 5;

        writer->initializeFile("Block Test File", num_dim,
        num_nodes, num_el, num_el_blk, num_ns, num_ss);

        std::string blockType = "EX_ELEM_BLOCK";
        int blockIdentifier = 10;
        std::string description = "quad";
        int numberOfEntries = 1;
        int numberOfNodesPerEntry = 4;
        int numberOfEdgesPerEntry = 4;
        int numberOfFacesPerEntry = 1;
        int numberOfAttributesPerEntry = 0;

        writer->defineBlock(blockType,blockIdentifier,description,numberOfEntries,
        numberOfNodesPerEntry,numberOfEdgesPerEntry,numberOfFacesPerEntry,
        numberOfAttributesPerEntry);


        int test_numberOfEntriesInBlock;
        int test_numberOfNodesPerEntry;
        int test_numberOfEdgesPerEntry;
        int test_numberOfFacesPerEntry;
        ex_get_block(writer->getExodusFileIdentifier(),EX_ELEM_BLOCK,blockIdentifier,NULL,&test_numberOfEntriesInBlock,&test_numberOfNodesPerEntry,&test_numberOfEdgesPerEntry,&test_numberOfFacesPerEntry,NULL);
        assert(test_numberOfEntriesInBlock == 1);
        assert(test_numberOfNodesPerEntry == 4);
        assert(test_numberOfEdgesPerEntry == 4);
        assert(test_numberOfFacesPerEntry == 1);


        writer->closeFile();
        std::remove("blockTest.exo");
}

void test_ExodusFileWriter_canSetBlockName(){
        std::remove("blockNameTest.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("blockNameTest.exo");
        int num_dim, num_nodes, num_el, num_el_blk, num_ns, num_ss;

        num_dim = 3; num_nodes = 33;
        num_el = 7; num_el_blk = 7;
        num_ns = 2; num_ss = 5;

        writer->initializeFile("Block Name Test File", num_dim,
        num_nodes, num_el, num_el_blk, num_ns, num_ss);

        std::string blockType = "EX_ELEM_BLOCK";
        int blockIdentifier = 10;
        std::string description = "quad";
        int numberOfEntries = 1;
        int numberOfNodesPerEntry = 4;
        int numberOfEdgesPerEntry = 0;
        int numberOfFacesPerEntry = 0;
        int numberOfAttributesPerEntry = 1;

        writer->defineBlock(blockType,blockIdentifier,description,numberOfEntries,
        numberOfNodesPerEntry,numberOfEdgesPerEntry,numberOfFacesPerEntry,
        numberOfAttributesPerEntry);

        writer->setBlockName(blockIdentifier,"Magic Block");

        writer->closeFile();
        std::remove("blockNameTest.exo");
}

void test_ExodusFileWriter_cantSetBlockNameWithoutBlockExisting(){
        std::remove("blockNameTestNonexistent.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("blockNameTestNonexistent.exo");
        int num_dim, num_nodes, num_el, num_el_blk, num_ns, num_ss;

        num_dim = 3; num_nodes = 33;
        num_el = 7; num_el_blk = 7;
        num_ns = 2; num_ss = 5;

        writer->initializeFile("Block Name Test File", num_dim,
        num_nodes, num_el, num_el_blk, num_ns, num_ss);


        try {
                int blockIdentifier = 372;
                writer->setBlockName(blockIdentifier,"Name of Nonexistent Block");
        } catch(std::exception& e){
                writer->closeFile();
                std::remove("blockNameTestNonexistent.exo");
                return;
        }
        writer->closeFile();
        std::remove("blockNameTestNonexistent.exo");
        assert(false && "The call in the try block should have generated an exception.");

}


void test_ExodusFileWriter_canSetBlockConnectivity(){
        std::remove("blockConnectivityTest.exo");
        int numberOfDimensions = 2;
        int numberOfNodes = 4;
        int numberOfElements = 1;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 1;
        int numberOfSideSets = 1;

        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();

        writer->openFile("blockConnectivityTest.exo");
        writer->initializeFile("Block Connectivity Test", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks,
        numberOfNodeSets, numberOfSideSets);

        std::string blockType = "EX_ELEM_BLOCK";
        int blockIdentifier = 932;
        std::string description = "quad";
        int numberOfEntries = 1;
        int numberOfNodesPerEntry = 4;
        int numberOfEdgesPerEntry = 0;
        int numberOfFacesPerEntry = 0;
        int numberOfAttributesPerEntry = 0;


        writer->defineBlock(blockType,blockIdentifier,description,numberOfEntries,
        numberOfNodesPerEntry,numberOfEdgesPerEntry,numberOfFacesPerEntry,
        numberOfAttributesPerEntry);

        std::list<int> nodeConnectivityArray{{1,2,3,4}};
        writer->setBlockConnectivity(blockType,blockIdentifier,&nodeConnectivityArray,NULL,NULL);
        writer->closeFile();
        std::remove("blockConnectivityTest.exo");
}

void test_ExodusFileWriter_setBlockConnectivityBlockMustExist(){
        std::remove("blockConnectivityBlockMustExistTest.exo");
        int numberOfDimensions = 2;
        int numberOfNodes = 4;
        int numberOfElements = 1;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 1;
        int numberOfSideSets = 1;

        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();

        writer->openFile("blockConnectivityBlockMustExistTest.exo");
        writer->initializeFile("Block Connectivity Block Must Exist Test", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks,
        numberOfNodeSets, numberOfSideSets);

        std::string blockType = "EX_ELEM_BLOCK";
        int blockIdentifier = 932;
        std::string description = "quad";
        int numberOfEntries = 1;
        int numberOfNodesPerEntry = 4;
        int numberOfEdgesPerEntry = 0;
        int numberOfFacesPerEntry = 0;
        int numberOfAttributesPerEntry = 0;

        try {
                std::list<int> nodeConnectivityArray{{1,2,3,4}};
                writer->setBlockConnectivity(blockType,blockIdentifier,&nodeConnectivityArray,NULL,NULL);
        } catch(std::exception& e){
                std::remove("blockConnectivityBlockMustExistTest.exo");
                return;
        }
        std::remove("blockConnectivityBlockMustExistTest.exo");
        assert(false && "The call in the try block should have generated an exception.");
}

void test_ExodusFileWriter_setBlockConnectivityCantReferenceExcessEntities(){
        std::remove("blockConnectivityExcessTest.exo");
        int numberOfDimensions = 2;
        int numberOfNodes = 4;
        int numberOfElements = 1;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 1;
        int numberOfSideSets = 1;

        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();

        writer->openFile("blockConnectivityExcessTest.exo");
        writer->initializeFile("Block Connectivity Excess Entities Test", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks,
        numberOfNodeSets, numberOfSideSets);

        std::string blockType = "EX_ELEM_BLOCK";
        int blockIdentifier = 932;
        std::string description = "quad";
        int numberOfEntries = 1;
        int numberOfNodesPerEntry = 4;
        int numberOfEdgesPerEntry = 0;
        int numberOfFacesPerEntry = 0;
        int numberOfAttributesPerEntry = 0;


        writer->defineBlock(blockType,blockIdentifier,description,numberOfEntries,
        numberOfNodesPerEntry,numberOfEdgesPerEntry,numberOfFacesPerEntry,
        numberOfAttributesPerEntry);

        try {
                //There are only 4 nodes in the block, but 10 entries in the connectivity list.
                std::list<int> nodeConnectivityArray{{1,2,3,4,5,6,7,8,9,10}};
                writer->setBlockConnectivity(blockType,blockIdentifier,&nodeConnectivityArray,NULL,NULL);
        } catch(std::exception& e){
                writer->closeFile();
                std::remove("blockConnectivityExcessTest.exo");
                return;
        }
        writer->closeFile();
        std::remove("blockConnectivityExcessTest.exo");
        assert(false && "The call in the try block should have generated an exception.");
}




void test_ExodusFileWriter_canSetNodeMap(){
        std::remove("nodeMapTest.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("nodeMapTest.exo");
        int numberOfDimensions = 3;
        int numberOfNodes = 8;
        int numberOfElements = 2;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 0;
        int numberOfSideSets = 0;

        writer->initializeFile("Node Map Test", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks, numberOfNodeSets, numberOfSideSets);

        std::list<int> nodeMap{10,20,30,40,50,60,70,80};
        writer->setNodeMap(&nodeMap);
        writer->closeFile();
        std::remove("nodeMapTest.exo");
}


void test_ExodusFileWriter_cantSetNodeMapBeforeOpeningFile(){
        std::remove("nodeMapTestNotOpen.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        int numberOfDimensions = 3;
        int numberOfNodes = 8;
        int numberOfElements = 2;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 0;
        int numberOfSideSets = 0;

        try {
                //Not Open!
                std::list<int> nodeMap{{10,20,30,40,50,60,70,80}};
                writer->setNodeMap(&nodeMap);
        } catch(std::exception& e){
                std::remove("nodeMapTestNotOpen.exo");
                return;
        }
        std::remove("nodeMapTestNotOpen.exo");
        assert(false && "The call in the try block should have generated an exception.");
}



void test_ExodusFileWriter_cantHaveEmptyNodeMap(){
        std::remove("nodeMapTestEmpty.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("nodeMapTestEmpty.exo");
        int numberOfDimensions = 3;
        int numberOfNodes = 8;
        int numberOfElements = 2;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 0;
        int numberOfSideSets = 0;

        writer->initializeFile("Node Map Test Empty", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks, numberOfNodeSets, numberOfSideSets);

        try {
                std::list<int> nodeMap;
                writer->setNodeMap(&nodeMap);
        } catch(std::exception& e){
                std::remove("nodeMapTestEmpty.exo");
                return;
        }
        std::remove("nodeMapTestEmpty.exo");
        assert(false && "The call in the try block should have generated an exception.");
}

/*
void test_ExodusFileWriter_cantHaveInvalidNodeMap(){
        std::remove("nodeMapTestInvalid.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("nodeMapTestInvalid.exo");
        int numberOfDimensions = 3;
        int numberOfNodes = 8;
        int numberOfElements = 2;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 0;
        int numberOfSideSets = 0;

        writer->initializeFile("Node Map Test Invalid", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks, numberOfNodeSets, numberOfSideSets);

        try {
                //Invalid! There are more node mappings than there are nodes (8 vs. 13)
                std::list<int> nodeMap{{10,20,30,40,50,60,70,80,90,100,110,120,130}};
                writer->setNodeMap(&nodeMap);
        } catch(std::exception& e){
                std::remove("nodeMapTestInvalid.exo");
                return;
        }
        std::remove("nodeMapTestInvalid.exo");
        assert(false && "The call in the try block should have generated an exception.");
}
*/

/*
void test_ExodusFileWriter_cantSetNodeMapMoreThanOnce(){
        std::remove("nodeMapTestRepeat.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("nodeMapTestRepeat.exo");
        int numberOfDimensions = 3;
        int numberOfNodes = 8;
        int numberOfElements = 2;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 0;
        int numberOfSideSets = 0;

        writer->initializeFile("Node Map Test Repeat", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks, numberOfNodeSets, numberOfSideSets);

        try {
                std::list<int> nodeMap{{10,20,30,40,50,60,70,80}};
                writer->setNodeMap(&nodeMap);
                writer->setNodeMap(&nodeMap);
        } catch(std::exception& e){
                std::remove("nodeMapTestRepeat.exo");
                return;
        }
        std::remove("nodeMapTestRepeat.exo");
        assert(false && "The call in the try block should have generated an exception.");
}
*/



void test_ExodusFileWriter_canSetElementMap(){
        std::remove("elementMapTest.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("elementMapTest.exo");
        int numberOfDimensions = 3;
        int numberOfNodes = 8;
        int numberOfElements = 2;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 0;
        int numberOfSideSets = 0;

        writer->initializeFile("Element Map Test", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks, numberOfNodeSets, numberOfSideSets);

        std::list<int> elementMap{{100,101}};
        writer->setElementMap(&elementMap);
        std::remove("elementMapTest.exo");
}

void test_ExodusFileWriter_cantSetElementMapBeforeOpeningFile(){
        std::remove("elementMapTestNotOpen.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        int numberOfDimensions = 3;
        int numberOfNodes = 8;
        int numberOfElements = 2;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 0;
        int numberOfSideSets = 0;

        try {
                //Not Open!
                std::list<int> elementMap{{10,20,30,40,50,60,70,80}};
                writer->setElementMap(&elementMap);
        } catch(std::exception& e){
                std::remove("elementMapTestNotOpen.exo");
                return;
        }
        std::remove("elementMapTestNotOpen.exo");
        assert(false && "The call in the try block should have generated an exception.");
}



void test_ExodusFileWriter_cantHaveEmptyElementMap(){
        std::remove("elementMapTestEmpty.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("elementMapTestEmpty.exo");
        int numberOfDimensions = 3;
        int numberOfNodes = 8;
        int numberOfElements = 2;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 0;
        int numberOfSideSets = 0;

        writer->initializeFile("Element Map Test Empty", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks, numberOfNodeSets, numberOfSideSets);

        try {
                std::list<int> elementMap;
                writer->setElementMap(&elementMap);
        } catch(std::exception& e){
                std::remove("elementMapTestEmpty.exo");
                return;
        }
        std::remove("elementMapTestEmpty.exo");
        assert(false && "The call in the try block should have generated an exception.");
}

void test_ExodusFileWriter_canDefineNodeSet(){
        std::remove("defineNodeSetTest.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("defineNodeSetTest.exo");
        int numberOfDimensions = 3;
        int numberOfNodes = 8;
        int numberOfElements = 2;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 2;
        int numberOfSideSets = 0;

        writer->initializeFile("Define Node Set Test", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks, numberOfNodeSets, numberOfSideSets);

        std::list<int> nodeSetA{{1,2,3,4}};
        int nodeSetAIdentifier = 80;
        std::list<int> nodeSetB{{5,6,7,8}};
        int nodeSetBIdentifier = 81;
        writer->defineNodeSet(nodeSetAIdentifier,"Node Set A",&nodeSetA);
        writer->defineNodeSet(nodeSetBIdentifier,"Node Set B",&nodeSetB);
        writer->closeFile();
        std::remove("defineNodeSetTest.exo");
}

void test_ExodusFileWriter_canDefineSideSet(){
        std::remove("defineSideSetTest.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("defineSideSetTest.exo");
        int numberOfDimensions = 3;
        int numberOfNodes = 8;
        int numberOfElements = 2;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 0;
        int numberOfSideSets = 2;

        writer->initializeFile("Define Side Set Test", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks, numberOfNodeSets, numberOfSideSets);

        std::list<int> elementSetA{{1,2,3,4}};
        std::list<int> sideSetA{{1,1,1,1}};
        int sideSetAIdentifier = 80;
        std::list<int> elementSetB{{5,6,7,8}};
        std::list<int> sideSetB{{1,1,1,1}};
        int sideSetBIdentifier = 81;
        writer->defineSideSet(sideSetAIdentifier,"Side Set A",&elementSetA,&sideSetA);
        writer->defineSideSet(sideSetBIdentifier,"Side Set B",&elementSetB,&sideSetB);
        writer->closeFile();
        std::remove("defineSideSetTest.exo");
}

void test_ExodusFileWriter_writeSquare(){
        std::remove("square.exo");
        int numberOfDimensions = 2; //This is a two-dimensional object.
        int numberOfNodes = 4; //The object has four nodes.
        int numberOfElements = 1; //These four nodes are grouped into one quad4 element.
        int numberOfElementBlocks = 1; //That element is grouped under a single element block.
        int numberOfNodeSets = 1; //There are no node sets.
        int numberOfSideSets = 1; //There are no side sets.

        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        ex_opts(EX_VERBOSE);

        writer->openFile("square.exo");
        writer->initializeFile("2D Square", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks,
        numberOfNodeSets, numberOfSideSets);


        std::list<double> xCoordinates{0.0,1.0,0.0,1.0};
        std::list<double> yCoordinates{1.0,0.0,0.0,1.0};

        std::list<std::string> coordinateNames = {"X", "Y"};
        writer->setCoordinatesOnAxis(&xCoordinates,0);
        writer->setCoordinatesOnAxis(&yCoordinates,1);
        writer->setCoordinateNames(coordinateNames);


        std::list<int> nodeIdentifierMap{10,20,30,40};
        writer->setNodeMap(&nodeIdentifierMap);
        std::list<int> elementIdentifierMap{100};
        writer->setElementMap(&elementIdentifierMap);

        std::string blockType = "EX_ELEM_BLOCK";
        int blockIdentifier = 101;
        std::string description = "quad";
        int numberOfEntries = 1;
        int numberOfNodesPerEntry = 4;
        int numberOfEdgesPerEntry = 0;
        int numberOfFacesPerEntry = 0;
        int numberOfAttributesPerEntry = 1;



        writer->defineBlock(blockType,blockIdentifier,description,numberOfEntries,
        numberOfNodesPerEntry,numberOfEdgesPerEntry,numberOfFacesPerEntry,
        numberOfAttributesPerEntry);

        std::list<int> nodeConnectivityArray{4,2,3,1};
        writer->setBlockConnectivity(blockType,blockIdentifier,&nodeConnectivityArray,NULL,NULL);


        std::string blockName = "The Square";
        writer->setBlockName(blockIdentifier,blockName);

        const char *attn_T[] = {"Peculiarity"};
        double      attr_T[] = {75.32};
        ex_put_attr(writer->getExodusFileIdentifier(), EX_ELEM_BLOCK, blockIdentifier, attr_T);

        std::list<int>nodeSetNodes{1, 2, 3, 4};
        int nodeSetId = 8000;
        writer->defineNodeSet(nodeSetId,"Quad_Node_Set",&nodeSetNodes);

        /* Write Side Set*/

        std::list<int>sideSetElements{1};
        std::list<int>sideSetSides{1};
        int sideSetId = 8001;
        writer->defineSideSet(sideSetId,"Quad_Element_Side_Set",&sideSetElements,&sideSetSides);

        writer->closeFile();

        /*
        (From page 14 of the documentation)
        As an example, suppose a database contains exactly one QUAD element with four nodes.
        The user desires the element ID to be 100 and the node IDs to be 10, 20, 30, and 40.


        The internal data structures representing the above model would be the following:
        - nodal coordinate array: (0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0)
        - connectivity array: (1, 2, 3, 4)
        - node number map: (10, 20, 30, 40)
        - element number map: (100)

        Internal (contiguously numbered) node and element IDs must be used for all data structures
        that contain node or element numbers (IDs), including node set node lists,
        side set element lists, and element connectivity.
        */


}



void test_ExodusFileWriter_writeTestExample2D(){
        /* Based on twod.c Exodus test. */

        std::remove("twod.exo");
        MeshingExodusFileWriter* writer = new MeshingExodusFileWriter();
        writer->openFile("twod.exo");
        int numberOfDimensions = 2;
        int numberOfNodes = 13;
        int numberOfElements = 20;
        int numberOfElementBlocks = 5;
        int numberOfNodeSets = 2;
        int numberOfSideSets = 2;

        ex_opts(EX_VERBOSE);

        writer->initializeFile("This is a 2D mesh example with tri, quad, beam, truss, circle", numberOfDimensions,
        numberOfNodes, numberOfElements, numberOfElementBlocks, numberOfNodeSets, numberOfSideSets);
        int exoid = writer->getExodusFileIdentifier();



        std::list<double> xCoordinates{0.0,-0.5,0.5,0.5,-0.5,-1.0,1.0,1.0,-1.0,-2.0,0.0,2.0,0.0};
        std::list<double> yCoordinates{0.0,-0.5,-0.5,0.5,0.5,-1.0,-1.0,1.0,1.0,0.0,-2.0,0.0,2.0};


        writer->setCoordinatesOnAxis(&xCoordinates, 0);
        writer->setCoordinatesOnAxis(&yCoordinates, 1);

        std::list<std::string> coordinateNames = {"xcoor", "ycoor"};
        writer->setCoordinateNames(coordinateNames);

        std::list<int> nodeMap{10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130};
        writer->setNodeMap(&nodeMap);

        std::list<int> elementMap{11,  21,  31,  41,  52,  62,  72,  82,  93,  103,
                              113, 123, 133, 143, 153, 163, 174, 184, 194, 204};
        writer->setElementMap(&elementMap);

        //const char *block_names[]        = {"Triangles", "Quadrilaterals", "", "Trusses", "Circles"};
        int         num_elem_in_block[]  = {4, 4, 4, 4, 4};
        int         num_nodes_per_elem[] = {3, 4, 2, 2, 1};
        int         ebids[]       = {100, 200, 300, 400, 500};

        writer->defineBlock("EX_ELEM_BLOCK", ebids[0], "triangle", num_elem_in_block[0],
                         num_nodes_per_elem[0], 0, 0, 0);

        writer->defineBlock("EX_ELEM_BLOCK", ebids[1], "quad", num_elem_in_block[1],
                         num_nodes_per_elem[1], 0, 0, 0);

        writer->defineBlock("EX_ELEM_BLOCK", ebids[2], "beam", num_elem_in_block[2],
                         num_nodes_per_elem[2], 0, 0, 3);

        writer->defineBlock("EX_ELEM_BLOCK", ebids[3], "truss", num_elem_in_block[3],
                         num_nodes_per_elem[3], 0, 0, 1);

        writer->defineBlock("EX_ELEM_BLOCK", ebids[4], "circle", num_elem_in_block[4],
                         num_nodes_per_elem[4], 0, 0, 2);

        writer->setBlockName(ebids[0],"Triangles");
        writer->setBlockName(ebids[1],"Quadrilaterals");
        writer->setBlockName(ebids[2],"");
        writer->setBlockName(ebids[3],"Trusses");
        writer->setBlockName(ebids[4],"Circles");

        std::list<int>conn_t{2, 3, 1, 3, 4, 1, 4, 5, 1, 5, 2, 1};
        std::list<int>conn_q{6, 7, 3, 2, 7, 8, 4, 3, 8, 9, 5, 4, 9, 6, 2, 5};
        std::list<int>conn_B{11, 7, 8, 13, 13, 9, 6, 11};
        std::list<int>conn_T{10, 6, 9, 10, 7, 12, 12, 8};
        std::list<int>conn_c{6, 7, 8, 9};

        writer->setBlockConnectivity("EX_ELEM_BLOCK",ebids[0],&conn_t,NULL,NULL);
        writer->setBlockConnectivity("EX_ELEM_BLOCK",ebids[1],&conn_q,NULL,NULL);
        writer->setBlockConnectivity("EX_ELEM_BLOCK",ebids[2],&conn_B,NULL,NULL);
        writer->setBlockConnectivity("EX_ELEM_BLOCK",ebids[3],&conn_T,NULL,NULL);
        writer->setBlockConnectivity("EX_ELEM_BLOCK",ebids[4],&conn_c,NULL,NULL);

          int nodeSetIDs[] = {20, 22};
          std::list<int> nod1{5, 4, 3, 2, 1};
          std::list<int> nod2{6, 7, 8, 9, 2, 3, 4, 5};

          writer->defineNodeSet(nodeSetIDs[0],"Triangle_Nodes",&nod1);
          writer->defineNodeSet(nodeSetIDs[1],"Quadrilateral_Nodes",&nod2);

          int sideSetIDs[] = {100,200};
          std::list<int> ss1el{1, 2, 3, 4};
          std::list<int> ss1si{1, 1, 1, 1};

          std::list<int> ss2el{5, 7, 6, 8};
          std::list<int> ss2si{1, 1, 1, 1};

          writer->defineSideSet(sideSetIDs[0],"A", &ss1el, &ss1si);
          writer->defineSideSet(sideSetIDs[1],"B", &ss2el, &ss2si);

          writer->closeFile();

}



void test_cube_error_reproduction(){

        /***

        The results should look like this:

                          .+------+
                         .' |    .'|  Cube #1, covered
                        +---+--+'  |  by face/element
                        |   |  |   |  blocks A/B
                        |  ,+--+---+
                        |.'    | .'
                        +------+'
          .+------+
         .' |    .'| Cube #2, covered
        +---+--+'  | by face/element
        |   |  |   | blocks A/B
        |  ,+--+---+
        |.'    | .'
        +------+'

        ***/

        int numberOfDimensions = 3;
        int numberOfNodes = 48;
        int numberOfElements = 2;
        int numberOfElementBlocks = 1;
        int numberOfNodeSets = 0;
        int numberOfSideSets = 0;
        int numberOfFaces = 12;
        int numberOfFaceBlocks = 1;

        float x[48], y[48], z[48];

        /* Cube #1 Nodes */
        x[0]=0.0; y[0]=0.0; z[0]=0.0;
        x[1]=0.0; y[1]=0.0; z[1]=1.0;
        x[2]=0.0; y[2]=1.0; z[2]=1.0;
        x[3]=0.0; y[3]=1.0; z[3]=0.0;
        x[4]=1.0; y[4]=1.0; z[4]=0.0;
        x[5]=1.0; y[5]=1.0; z[5]=1.0;
        x[6]=1.0; y[6]=0.0; z[6]=1.0;
        x[7]=1.0; y[7]=0.0; z[7]=0.0;
        x[8]=0.0; y[8]=0.0; z[8]=0.0;
        x[9]=1.0; y[9]=0.0; z[9]=0.0;
        x[10]=1.0; y[10]=0.0; z[10]=1.0;
        x[11]=0.0; y[11]=0.0; z[11]=1.0;
        x[12]=0.0; y[12]=0.0; z[12]=1.0;
        x[13]=1.0; y[13]=0.0; z[13]=1.0;
        x[14]=1.0; y[14]=1.0; z[14]=1.0;
        x[15]=0.0; y[15]=1.0; z[15]=1.0;
        x[16]=0.0; y[16]=1.0; z[16]=1.0;
        x[17]=1.0; y[17]=1.0; z[17]=1.0;
        x[18]=1.0; y[18]=1.0; z[18]=0.0;
        x[19]=0.0; y[19]=1.0; z[19]=0.0;
        x[20]=0.0; y[20]=1.0; z[20]=0.0;
        x[21]=1.0; y[21]=1.0; z[21]=0.0;
        x[22]=1.0; y[22]=0.0; z[22]=0.0;
        x[23]=0.0; y[23]=0.0; z[23]=0.0;

        /* Cube #2 Nodes */
        x[24]=5.0; y[24]=5.0; z[24]=5.0;
        x[25]=5.0; y[25]=5.0; z[25]=6.0;
        x[26]=5.0; y[26]=6.0; z[26]=6.0;
        x[27]=5.0; y[27]=6.0; z[27]=5.0;
        x[28]=6.0; y[28]=6.0; z[28]=5.0;
        x[29]=6.0; y[29]=6.0; z[29]=6.0;
        x[30]=6.0; y[30]=5.0; z[30]=6.0;
        x[31]=6.0; y[31]=5.0; z[31]=5.0;
        x[32]=5.0; y[32]=5.0; z[32]=5.0;
        x[33]=6.0; y[33]=5.0; z[33]=5.0;
        x[34]=6.0; y[34]=5.0; z[34]=6.0;
        x[35]=5.0; y[35]=5.0; z[35]=6.0;
        x[36]=5.0; y[36]=5.0; z[36]=6.0;
        x[37]=6.0; y[37]=5.0; z[37]=6.0;
        x[38]=6.0; y[38]=6.0; z[38]=6.0;
        x[39]=5.0; y[39]=6.0; z[39]=6.0;
        x[40]=5.0; y[40]=6.0; z[40]=6.0;
        x[41]=6.0; y[41]=6.0; z[41]=6.0;
        x[42]=6.0; y[42]=6.0; z[42]=5.0;
        x[43]=5.0; y[43]=6.0; z[43]=5.0;
        x[44]=5.0; y[44]=6.0; z[44]=5.0;
        x[45]=6.0; y[45]=6.0; z[45]=5.0;
        x[46]=6.0; y[46]=5.0; z[46]=5.0;
        x[47]=5.0; y[47]=5.0; z[47]=5.0;
}

void test_nfaced(){
        int exoid, num_dim, num_nodes, num_elem, num_elem_blk;
        int num_elem_in_block[10], num_total_nodes_per_blk[10];
        int num_face_in_block[10], num_total_faces_per_blk[10];
        int num_node_sets, error;
        int i, j, *connect;
        int bids, nnpe[20];
        int num_qa_rec, num_info;
        int CPU_word_size, IO_word_size;

        float x[100], y[100], z[100];
        char *block_names[10];

        std::string title = "This is a test";
        ex_opts(EX_VERBOSE | EX_ABORT);

        /* Specify compute and i/o word size */

        CPU_word_size = 0; /* sizeof(float) */
        IO_word_size  = 4; /* (4 bytes) */

        /* create EXODUS II file */

        exoid = ex_create("test-nfaced.exo", /* filename path */
                          EX_CLOBBER,        /* create mode */
                          &CPU_word_size,    /* CPU float word size in bytes */
                          &IO_word_size);    /* I/O float word size in bytes */
        printf("after ex_create for test.exo, exoid = %d\n", exoid);
        printf(" cpu word size: %d io word size: %d\n", CPU_word_size, IO_word_size);

        /* initialize file with parameters */
        {
          ex_init_params par;

          num_dim       = 3;
          num_nodes     = 14;
          num_elem      = 3;
          num_elem_blk  = 1;
          num_node_sets = 0;

          ex_copy_string(par.title, title.c_str(), MAX_LINE_LENGTH + 1);
          par.num_dim       = num_dim;
          par.num_nodes     = num_nodes;
          par.num_edge      = 0;
          par.num_edge_blk  = 0;
          par.num_face      = 15;
          par.num_face_blk  = 1;
          par.num_elem      = num_elem;
          par.num_elem_blk  = num_elem_blk;
          par.num_node_sets = num_node_sets;
          par.num_edge_sets = 0;
          par.num_face_sets = 0;
          par.num_side_sets = 0;
          par.num_elem_sets = 0;
          par.num_node_maps = 0;
          par.num_edge_maps = 0;
          par.num_face_maps = 0;
          par.num_elem_maps = 0;

          error = ex_put_init_ext(exoid, &par);

          printf("after ex_put_init_ext, error = %d\n", error);

          if (error) {
            ex_close(exoid);
            exit(-1);
          }
        }

        /* write nodal coordinates values and names to database */
        x[0]  = 0.00000e+00;
        y[0]  = 0.00000e+00;
        z[0]  = 0.00000e+00;
        x[1]  = 2.00000e+00;
        y[1]  = 0.00000e+00;
        z[1]  = 0.00000e+00;
        x[2]  = 0.00000e+00;
        y[2]  = 2.00000e+00;
        z[2]  = 0.00000e+00;
        x[3]  = 2.00000e+00;
        y[3]  = 2.00000e+00;
        z[3]  = 0.00000e+00;
        x[4]  = 0.00000e+00;
        y[4]  = 0.00000e+00;
        z[4]  = 2.00000e+00;
        x[5]  = 2.00000e+00;
        y[5]  = 0.00000e+00;
        z[5]  = 2.00000e+00;
        x[6]  = 0.00000e+00;
        y[6]  = 2.00000e+00;
        z[6]  = 2.00000e+00;
        x[7]  = 2.00000e+00;
        y[7]  = 2.00000e+00;
        z[7]  = 2.00000e+00;
        x[8]  = 0.00000e+00;
        y[8]  = 3.50000e+00;
        z[8]  = 1.00000e+00;
        x[9]  = 2.00000e+00;
        y[9]  = 3.50000e+00;
        z[9]  = 1.00000e+00;
        x[10] = 0.00000e+00;
        y[10] = 3.00000e+00;
        z[10] = 1.50000e+00;
        x[11] = 2.00000e+00;
        y[11] = 3.00000e+00;
        z[11] = 1.50000e+00;
        x[12] = 0.00000e+00;
        y[12] = 3.00000e+00;
        z[12] = 0.50000e+00;
        x[13] = 2.00000e+00;
        y[13] = 3.00000e+00;
        z[13] = 0.50000e+00;

        error = ex_put_coord(exoid, x, y, z);
        printf("after ex_put_coord, error = %d\n", error);

        if (error) {
          ex_close(exoid);
          exit(-1);
        }

        const char* coord_names[3] = { "x", "y", "z" };
        //coord_names[0] = "x";
        //coord_names[1] = "y";
        //coord_names[2] = "z";

        error = ex_put_coord_names(exoid, (char**)coord_names);
        printf("after ex_put_coord_names, error = %d\n", error);

        if (error) {
          ex_close(exoid);
          exit(-1);
        }

        /* Write the face block parameters */
        block_names[0]             = "face_block_1";
        num_face_in_block[0]       = 15;
        num_total_nodes_per_blk[0] = 58;
        bids                       = 10;

        error = ex_put_block(exoid, EX_FACE_BLOCK, bids, "nsided", num_face_in_block[0],
                             num_total_nodes_per_blk[0], 0, 0, 0);
        printf("after ex_put_block, error = %d\n", error);

        if (error) {
          ex_close(exoid);
          exit(-1);
        }

        /* Write face block names */
        error = ex_put_names(exoid, EX_FACE_BLOCK, block_names);
        printf("after ex_put_names, error = %d\n", error);

        if (error) {
          ex_close(exoid);
          exit(-1);
        }

        /* write face connectivity */

        connect = (int *)calloc(num_total_nodes_per_blk[0], sizeof(int));

        i = 0;
        j = 0;


        connect[i++] = 5;
        connect[i++] = 6;
        connect[i++] = 8; // connectivity of face 1 of element 1
        nnpe[j++]    = 3;

        connect[i++] = 2;
        connect[i++] = 1;
        connect[i++] = 4; // face 2 of element 1
        nnpe[j++]    = 3;



        connect[i++] = 6;
        connect[i++] = 2;
        connect[i++] = 4;
        connect[i++] = 8; // face 3 of element 1
        nnpe[j++]    = 4;



        connect[i++] = 8;
        connect[i++] = 4;
        connect[i++] = 1;
        connect[i++] = 5; /* face 4 of element 1 */
        nnpe[j++]    = 4;

        connect[i++] = 1;
        connect[i++] = 2;
        connect[i++] = 6;
        connect[i++] = 5; /*  face 5 of element 1 */
        nnpe[j++]    = 4;

        connect[i++] = 5;
        connect[i++] = 8;
        connect[i++] = 7; /* connectivity of face 1 of element 2 */
        nnpe[j++]    = 3;

        connect[i++] = 1;
        connect[i++] = 3;
        connect[i++] = 4; /*  face 2 of element 2 */
        nnpe[j++]    = 3;

        connect[i++] = 7;
        connect[i++] = 8;
        connect[i++] = 4;
        connect[i++] = 3; /*  face 3 of element 2 */
        nnpe[j++]    = 4;

        connect[i++] = 7;
        connect[i++] = 3;
        connect[i++] = 1;
        connect[i++] = 5; /*  face 4 of element 2 */
        nnpe[j++]    = 4;

        connect[i++] = 8;
        connect[i++] = 4;
        connect[i++] = 14;
        connect[i++] = 10;
        connect[i++] = 12; /* connectivity of face 1 of element 3 */
        nnpe[j++]    = 5;

        connect[i++] = 7;
        connect[i++] = 11;
        connect[i++] = 9;
        connect[i++] = 13;
        connect[i++] = 3; /*  face 2 of element 3 */
        nnpe[j++]    = 5;

        connect[i++] = 7;
        connect[i++] = 8;
        connect[i++] = 12;
        connect[i++] = 11; /* face 3 of element 3 */
        nnpe[j++]    = 4;

        connect[i++] = 11;
        connect[i++] = 12;
        connect[i++] = 10;
        connect[i++] = 9; /* face 4 of element 3 */
        nnpe[j++]    = 4;

        connect[i++] = 9;
        connect[i++] = 10;
        connect[i++] = 14;
        connect[i++] = 13; /*  face 5 of element 3 */
        nnpe[j++]    = 4;

        connect[i++] = 13;
        connect[i++] = 14;
        connect[i++] = 4;
        connect[i++] = 3; /* face 6 of element 3 */
        nnpe[j++]    = 4;

        assert(i == num_total_nodes_per_blk[0]);
        assert(j == num_face_in_block[0]);

        error = ex_put_conn(exoid, EX_FACE_BLOCK, bids, connect, NULL, NULL);
        printf("after ex_put_conn, error = %d\n", error);

        if (error) {
          ex_close(exoid);
          exit(-1);
        }

        free(connect);
        connect = NULL;

        error = ex_put_entity_count_per_polyhedra(exoid, EX_FACE_BLOCK, bids, nnpe);
        printf("after ex_put_entity_count_per_polyhedra, error = %d\n", error);

        if (error) {
          ex_close(exoid);
          exit(-1);
        }

        /* write element block parameters */
        block_names[0] = "nfaced_1";

        num_elem_in_block[0]       = 3;
        num_total_faces_per_blk[0] = 5 + 5 + 7;

        bids = 10;

        error = ex_put_block(exoid, EX_ELEM_BLOCK, bids, "nfaced", num_elem_in_block[0], 0, 0,
                             num_total_faces_per_blk[0], 0);
        printf("after ex_put_block, error = %d\n", error);

        if (error) {
          ex_close(exoid);
          exit(-1);
        }

        /* Write element block names */
        error = ex_put_names(exoid, EX_ELEM_BLOCK, block_names);
        printf("after ex_put_names, error = %d\n", error);

        if (error) {
          ex_close(exoid);
          exit(-1);
        }

        /* write element-face connectivity */
        connect = (int *)calloc(num_total_faces_per_blk[0], sizeof(int));

        i            = 0;
        j            = 0;
        connect[i++] = 1;
        connect[i++] = 2;
        connect[i++] = 3;
        connect[i++] = 4;
        connect[i++] = 5;
        nnpe[j++]    = 5; /* Number of faces per element 1 */

        connect[i++] = 4;
        connect[i++] = 6;
        connect[i++] = 7;
        connect[i++] = 8;
        connect[i++] = 9;
        nnpe[j++]    = 5; /* Number of faces per element 2 */

        connect[i++] = 8;
        connect[i++] = 10;
        connect[i++] = 11;
        connect[i++] = 12;
        connect[i++] = 13;
        connect[i++] = 14;
        connect[i++] = 15;
        nnpe[j++]    = 7; /* Number of faces per element 3 */

        assert(i == num_total_faces_per_blk[0]);
        assert(j == num_elem_in_block[0]);

        error = ex_put_conn(exoid, EX_ELEM_BLOCK, bids, NULL, NULL, connect);
        printf("after ex_put_conn, error = %d\n", error);

        if (error) {
          ex_close(exoid);
          exit(-1);
        }

        free(connect);

        error = ex_put_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK, bids, nnpe);
        printf("after ex_put_entity_count_per_polyhedra, error = %d\n", error);

        if (error) {
          ex_close(exoid);
          exit(-1);
        }

        /* write QA records; test empty and just blank-filled records */
        num_qa_rec = 2;

        const char* qa_record[2][4] = {
          "TESTWT-NFACED",
          "testwt-nfaced",
          "2010/02/15",
          "06:35:15",
          "",
          "                            ",
          "",
          "                        ",
        };
        //qa_record[0][0] = "TESTWT-NFACED";
        //qa_record[0][1] = "testwt-nfaced";
        //qa_record[0][2] = "2010/02/15";
        //qa_record[0][3] = "06:35:15";
        //qa_record[1][0] = "";
        //qa_record[1][1] = "                            ";
        //qa_record[1][2] = "";
        //qa_record[1][3] = "                        ";

        error = ex_put_qa(exoid, num_qa_rec, (char*(*)[4])qa_record);
        printf("after ex_put_qa, error = %d\n", error);

        if (error) {
          ex_close(exoid);
          exit(-1);
        }

        /* write information records; test empty and just blank-filled records */
        num_info = 3;

        const char *info[3] = {
                "This is the first information record.",
                "",
                "                                     "
        };
        //info[0] = "This is the first information record.";
        //info[1] = "";
        //info[2] = "                                     ";

        error = ex_put_info(exoid, num_info, (char**)info);
        printf("after ex_put_info, error = %d\n", error);

        if (error) {
          ex_close(exoid);
          exit(-1);
        }

        /* close the EXODUS files
         */
        error = ex_close(exoid);
        printf("after ex_close, error = %d\n", error);
        if (error) {
          ex_close(exoid);
          exit(-1);
        }

}



#endif //USE_EXODUS

int main(int argc, char* argv[])
{
        #ifdef USE_EXODUS
        //Test MeshingExodusFileWriter

        //Instantiation and File Creation Tests
        test_ExodusFileWriter_isDirectlyConstructible();
        test_ExodusFileWriter_isDirectlyDestructible();
        test_ExodusFileWriter_DoesNotDetectNonExistentFile();
        test_ExodusFileWriter_canOpenNewFile();
        test_ExodusFileWriter_canCloseFile();
        test_ExodusFileWriter_closingNonexistentFileCausesException();
        test_ExodusFileWriter_canOpenAndCloseArbitrarilyManyFiles();

        //File Setup Tests
        test_ExodusFileWriter_canInitializeFile();
        test_ExodusFileWriter_cantCallInitializeFileMoreThanOnce();
        test_ExodusFileWriter_cantCallInitializeFileWithNoOpenFile();
        test_ExodusFileWriter_canSetCoordinateNames();
        test_ExodusFileWriter_mustInitializeFileBeforeSettingCoordinateNames();

        //Entity Tests
        test_ExodusFileWriter_canCreateBlock();
        test_ExodusFileWriter_canSetBlockName();
        test_ExodusFileWriter_cantSetBlockNameWithoutBlockExisting();

        test_ExodusFileWriter_canSetBlockConnectivity();
        //test_ExodusFileWriter_setBlockConnectivityBlockMustExist();
        //test_ExodusFileWriter_setBlockConnectivityCantReferenceExcessEntities();

        test_ExodusFileWriter_canSetNodeMap();
        test_ExodusFileWriter_cantSetNodeMapBeforeOpeningFile();
        test_ExodusFileWriter_cantHaveEmptyNodeMap();

        test_ExodusFileWriter_canSetElementMap();
        test_ExodusFileWriter_cantSetElementMapBeforeOpeningFile();
        test_ExodusFileWriter_cantHaveEmptyElementMap();

        test_ExodusFileWriter_canDefineNodeSet();
        test_ExodusFileWriter_canDefineSideSet();


        test_ExodusFileWriter_writeSquare();
        test_ExodusFileWriter_writeTestExample2D();



        test_nfaced();

        #endif //USE_EXODUS

        //Either USE_EXODUS or NO_EXODUS must be defined or this test will fail.
        #ifndef USE_EXODUS
                #ifndef NO_EXODUS
                        throw std::runtime_error("Either NO_EXODUS or USE_EXODUS must be defined!");
                #endif
        #endif

        return 0;
}


