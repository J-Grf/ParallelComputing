#include "tri.h"
#include <iostream>   // std::cout
#include <string>     // std::string, std::stoi
#include <sstream>
#include <math.h>       /* atan */

#define PI 3.14159265
/***************************************************************************************************
void preProcessor::prepareMesh()
****************************************************************************************************
Does everything related to mesh generation and transfers between processors.
***************************************************************************************************/
void triMesh::prepareMesh(inputSettings* settings)
{
    cout << endl << "====================== MESH ======================" << endl;

    readMeshFiles(settings);

    return;
}

/***************************************************************************************************
void preProcessor::readMeshFiles()
****************************************************************************************************
File read procedure :
1- Name of the file to be opened is retrieved from the inputSetting obj.
2- File is opened in appropriate format, this is ascii format for minf or
   other text files and binary format for binary mesh files.
3- Read operation for minf file is straight forward. Binary files are read
   as size of a double or int and stored in readStream. Then swapbytes
   function is called to swap the bytes for the correct endianness.
4- Finally obtained data is deep-copied to the mesh data structure.
***************************************************************************************************/
void triMesh::readMeshFiles(inputSettings* settings)
{
    ifstream file;        // file name object for serial file read
    string   dummy;       // dummy string to hold file names etc.
    char*    readStream;  // temperory var used for strings read from files
    double   dummyDouble; // temperory var used for double values read from files
    int npes = std::stoi(settings->getnpesStr()); // number of processors from settings
    const string partitionType = settings->getPartitionType();
    stringstream stst;

    /***********************************************************************************************/
    //COMPLETE TASK 1 HERE
    dummy = settings->getMinfFile();
    const string minf = dummy;
    cout << "Reading " << minf << endl;
    file.open(minf.c_str(), ios::in);
    
    if (file.is_open()==false){
        cout << "Unable to open input file: " << minf << " ! Aborting... " << endl;
        exit(0);
    }
    
    file >> dummy >> nn;
    file >> dummy >> ne;
    cout << "number of nodes: " << nn << endl
         << "number of elements: " << ne << endl;
   
    file.close(); 
    /***********************************************************************************************/

    //Allocation of memory for the mesh data structure
    xyz = new double [nn*nsd];
    node = new triNode[nn];
    elem = new triElement[ne];
    pID  = new int [nn];

    /**********************************************************************************************/
    // COMPLETE TASK 2 HERE
    dummy = settings->getMxyzFile();
    const string mxyz = dummy;
    cout << "Reading " << mxyz << endl;
    file.open(mxyz.c_str(), ios::in | ios::binary | ios::ate);
    
    if (file.is_open()==false){
        cout << "Unable to open input file: " << mxyz << " ! Aborting... " << endl;
        exit(0);
    }

    readStream = new char [nsd * sizeof(double)];
    file.seekg(0, ios::beg); //set reading position to beginning
    for(int i = 0; i < nn; i++){
        file.read(readStream, nsd * sizeof(double));
        swapBytes(readStream, nsd, sizeof(double));
        for(int j = 0; j < nsd; j++){
            dummyDouble = *((double*)readStream + j); 
            //readStream is array (pointer type to char), it is casted to double array
            xyz[i * nsd + j] = dummyDouble;
        }
        node[i].setX(xyz[i * nsd]);
        node[i].setY(xyz[i * nsd + 1]);
    }
    file.close();
    free(readStream);
    /**********************************************************************************************/

    /**********************************************************************************************/
    //COMPLETE TASK 3 HERE
    dummy = settings->getMienFile();
    const string mien = dummy;
    cout << "Reading " << mien << endl;
    file.open(mien.c_str(), ios::in | ios::binary | ios::ate);
    
    if (file.is_open()==false){
        cout << "Unable to open input file: " << mien << " ! Aborting... " << endl;
        exit(0);
    }

    readStream = new char [nen * sizeof(int)];
    file.seekg(0, ios::beg); //set reading position to beginning

    for(int i = 0; i < ne; i++){
        file.read(readStream, nen * sizeof(int));
        swapBytes(readStream, nen, sizeof(int));
        for(int j = 0; j < nen; j++){
            elem[i].setConn(j, *((int*)readStream + j) - 1);
        }
    }

    file.close();
    free(readStream);
    /**********************************************************************************************/


    /***********************************************************************************************
     * BONUS TASK: COMPLETE BT HERE
    ***********************************************************************************************/
    cout << "Number of processors for parallelisation of mesh generation: " << npes << endl;
    register int i; //iteration index for all the nodes
    register int ipe; // iteration index for all the pes
    
    // Round Robin Partitioning:
    int dist = 0;
    double PhiP;
    if (partitionType == "angle") 
        PhiP = 2 * PI / npes;
        cout << "PhiP: " << PhiP << endl;
    if (partitionType != "none") 
        cout << "partitionType: " << partitionType << endl;
    for(i = 0; i < nn; i++){
        // Round Robin Partitioning:
        if(partitionType == "robin")
            dist = i % npes;
        // Block Partitioning based on node number:
        else if (partitionType == "node")
            dist = i * npes / nn;
        // Block Partitioning based on angle:
        else if (partitionType == "angle"){
            int angle;
            angle = round(atan (node[i].getY()/node[i].getX()) * 180 / PI);
            if (node[i].getY() < 0){
                angle = (angle + 180);
            }else if(node[i].getX() < 0){
                angle = (angle + 360);
            }
            dist = angle * npes / 360;

            /* double angle;
            if(node[i].getY() > 0 && node[i].getX() > 0)
                angle = atan(abs(node[i].getY()) / abs(node[i].getX()));
            else if(node[i].getY() > 0 || node[i].getX() < 0 )
                angle = PI / 2 - atan(abs(node[i].getY()) / abs(node[i].getX())) + PI / 2; 
            else if(node[i].getY() < 0 || node[i].getX() < 0 ) 
                angle = atan(abs(node[i].getY()) / abs(node[i].getX())) + PI;
            else if(node[i].getY() < 0 || node[i].getX() > 0 ) 
                angle = PI / 2 - atan(abs(node[i].getY()) / abs(node[i].getX())) + 3/2 * PI;

            cout << "angle: " << angle << endl;
            for(int j = 0; j < npes; j++){
                if((j * PhiP) < angle && angle < ((j+1) * PhiP)){
                    //cout << "bounds: " << (j * PhiP) << " and " << ((j+1) * PhiP) << endl;
                    dist = j;
                }
            } */
        }else
            break;
        
        pID[i] = dist; // assign node to a processor
        node[i].setpeID((dist));
    }
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