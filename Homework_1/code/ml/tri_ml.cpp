
#include "tri.h"
#include <iostream>   // std::cout
#include <string>     // std::string, std::stoi
#include <sstream>
#include <math.h>       /* atan */

#define PI 3.14159265

/***************************************************************************************************
 * void preProcessor::prepareMesh()
 ****************************************************************************************************
 * Does everythin' realted to mesh generation and transfers between processors.
 ***************************************************************************************************/
void triMesh::prepareMesh(inputSettings* settings)
{
    cout << endl << "====================== MESH ======================" << endl;

    readMeshFiles(settings);

    return;
}

/***************************************************************************************************
 * void preProcessor::readMeshFiles()
 ****************************************************************************************************
 * File read procedure :
 * 1- Name of the file to be opened is retrieved from the inputSetting obj.
 * 2- File is opened in appropriate format, this is ascii format for minf or
 * other text files and binary format for binary mesh files.
 * 3- Read operation for minf file is straight forward. Binary files are read
 * as size of a double or int and stored in readStream. Then swapbytes
 * function is called to swap the bytes for the correct endianness.
 * 4- Finally obtained data is deep-copied to the mesh data structure.
 ***************************************************************************************************/
void triMesh::readMeshFiles(inputSettings* settings)
{
    ifstream    file;           // file name object for serial file read
    string      dummy;          // dummy string to hold file names etc.
    char*       readStream;     // temperory var used for strings read from files
    double      dummyDouble;    // temperory var used for double values read from files

    /***********************************************************************************************
     * READ THE MINF FILE
     * This file should hold the number of elements, nodes, space dimensions, element nodes and
     * element faces.
     ***********************************************************************************************/
    dummy = settings->getMinfFile();
    file.open(dummy.c_str(), ios::in);
    if (file.is_open()==false)
    {
        cout << "Unable to open minf file: ! Aborting... " << endl;
        exit(0);
    }
    file >> dummy >> nn;
    file >> dummy >> ne;
    cout << "> Number of mesh elements : " << ne << endl;
    cout << "> Number of nodes : " << nn << endl;
    cout << "> File read complete: minf" << endl;
    file.close();

    //Allocation of memory for the mesh data structure
    xyz = new double [nn*nsd];
    node = new triNode[nn];
    elem = new triElement[ne];
    pID  = new int [nn];

    /***********************************************************************************************
     * READ THE MXYZ FILE
     * This file contains the node coordinates
     ***********************************************************************************************/
    dummy = settings->getMxyzFile();
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }

    readStream = new char [nsd*sizeof(double)];
    file.seekg (0, ios::beg);

    for(int i=0; i<nn; i++)
    {
        file.read(readStream, nsd*sizeof(double));
        swapBytes(readStream, nsd, sizeof(double));
        node[i].setX(*((double*)readStream)); 
        node[i].setY(*((double*)readStream+1));
        xyz[i*nsd+xsd] = *((double*)readStream);
        xyz[i*nsd+ysd] = *((double*)readStream+1);
    }
    cout << "> File read complete: " << dummy << endl;
    file.close();
    free(readStream);

    /***********************************************************************************************
     * READ THE MIEN FILE
     * This file contains the element connectivity
     ***********************************************************************************************/
    dummy = settings->getMienFile();
    file.open(dummy.c_str(), ios::in|ios::binary|ios::ate);
    if (file.is_open()==false)
    {
        cout << "Unable to open file : " << dummy << endl;
        exit(0);
    }
    int connValue;
    readStream = new char [nen*sizeof(int)];
    file.seekg (0, ios::beg);

    for(int i=0; i<ne; i++)
    {
        file.read (readStream, nen*sizeof(int));
        swapBytes(readStream, nen, sizeof(int));
        for(int j=0; j<nen; j++)
        {
            connValue = *((int*)readStream+j);
            elem[i].setConn(j, connValue-1);
        }
    }

    cout << "> File read complete: " << dummy << endl;
    file.close();
    free(readStream);

    return;
}

void triMesh::swapBytes (char *array, int nelem, int elsize)
{
    register int sizet, sizem, i, j;
    char *bytea, *byteb;
    sizet = elsize;
    sizem = sizet - 1;
    bytea = new char [sizet];
    byteb = new char [sizet];
    for (i = 0; i < nelem; i++)
    {
        memcpy((void *)bytea, (void *)(array+i*sizet), sizet);
        for (j = 0; j < sizet; j++)
            byteb[j] = bytea[sizem - j];
        memcpy((void *)(array+i*sizet), (void *)byteb, sizet);
    }
    free(bytea);
    free(byteb);

    return;
}