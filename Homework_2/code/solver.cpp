#include "solver.h"
#include "settings.h"
/***************************************************************************************************
solverControl
****************************************************************************************************
Flow control to prepare and solve the system.
***************************************************************************************************/
void femSolver::solverControl(inputSettings* argSettings, triMesh* argMesh)
{
    cout << endl << "=================== SOLUTION =====================" << endl;

    mesh = argMesh;
    settings = argSettings;
    int ne = mesh->getNe();

    //loop over all elements
    for(int iec=0; iec<ne; iec++)
    {
        calculateJacobian(iec);
        calculateElementMatrices(iec);
    }
    applyDrichletBC();
    accumulateMass();
    explicitSolver();

    return;
}

/***************************************************************************************************
calculateJacobian
****************************************************************************************************
Compute and store the jacobian for each element.
***************************************************************************************************/
void femSolver::calculateJacobian(const int e)
{
    int myNode;     // node number for the current node.
    double x[nen];  // x values for all the nodes of an element.
    double y[nen];  // y values for all the nodes of an element.

    triMasterElement* ME = mesh->getME(0);  // for easy access to the master element. pointer to first ME GQ point
    double * xyz = mesh->getXyz();

    double J[2][2];     // Jacobian for the current element.
    double detJ;        // Jacobian determinant for the current element.
    double invJ[2][2];  // inverse of Jacobian for the current element.
    double dSdX[3];     // dSdx on a GQ point.
    double dSdY[3];     // dSdy on a GQ point.

    // collect element node coordinates in x[3] and y[3] matrices.
    // thetrahedric elements -> 3 nodes per element
    for (int i=0; i<nen; i++)
    {
        myNode =  mesh->getElem(e)->getConn(i);
        x[i] = xyz[myNode*nsd+xsd];
        y[i] = xyz[myNode*nsd+ysd];
    }

    // for all GQ points detJ, dSDx[3] and dSdY[3] are determined.
    for (int p=0; p<nGQP; p++)
    {
        // Calculate Jacobian
        J[0][0] = ME[p].getDSdKsi(0)*x[0] + ME[p].getDSdKsi(1)*x[1] + ME[p].getDSdKsi(2)*x[2];
        J[0][1] = ME[p].getDSdKsi(0)*y[0] + ME[p].getDSdKsi(1)*y[1] + ME[p].getDSdKsi(2)*y[2];
        J[1][0] = ME[p].getDSdEta(0)*x[0] + ME[p].getDSdEta(1)*x[1] + ME[p].getDSdEta(2)*x[2];
        J[1][1] = ME[p].getDSdEta(0)*y[0] + ME[p].getDSdEta(1)*y[1] + ME[p].getDSdEta(2)*y[2];

        //Calculate determinant of Jacobian and store in mesh.
        detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        mesh->getElem(e)->setDetJ(p, detJ);

        // Calculate inverse of Jacobian.
        invJ[0][0] =  J[1][1]/detJ;
        invJ[0][1] = -J[0][1]/detJ;
        invJ[1][0] = -J[1][0]/detJ;
        invJ[1][1] =  J[0][0]/detJ;

        // Calculate dSdx and dSdy and store in mesh.
        dSdX[0] = invJ[0][0]*ME[p].getDSdKsi(0) + invJ[0][1]*ME[p].getDSdEta(0);
        dSdX[1] = invJ[0][0]*ME[p].getDSdKsi(1) + invJ[0][1]*ME[p].getDSdEta(1);
        dSdX[2] = invJ[0][0]*ME[p].getDSdKsi(2) + invJ[0][1]*ME[p].getDSdEta(2);
        dSdY[0] = invJ[1][0]*ME[p].getDSdKsi(0) + invJ[1][1]*ME[p].getDSdEta(0);
        dSdY[1] = invJ[1][0]*ME[p].getDSdKsi(1) + invJ[1][1]*ME[p].getDSdEta(1);
        dSdY[2] = invJ[1][0]*ME[p].getDSdKsi(2) + invJ[1][1]*ME[p].getDSdEta(2);

        mesh->getElem(e)->setDSdX(p, 0, dSdX[0]);
        mesh->getElem(e)->setDSdX(p, 1, dSdX[1]);
        mesh->getElem(e)->setDSdX(p, 2, dSdX[2]);
        mesh->getElem(e)->setDSdY(p, 0, dSdY[0]);
        mesh->getElem(e)->setDSdY(p, 1, dSdY[1]);
        mesh->getElem(e)->setDSdY(p, 2, dSdY[2]);
    }

    return;
}

/***************************************************************************************************
void femSolver::calculateElementMatrices(const int e)
****************************************************************************************************
Compute the K, M and F matrices. Then accumulate the total mass into the node structure.
***************************************************************************************************/
void femSolver::calculateElementMatrices(const int e)
{
    int node;
    int D = settings->getD(); //diffusion coefficient / thermal diffusivity
    double f = settings->getSource(); //sourceTerm
    double *xyz = mesh->getXyz();

    double totalM = 0.0;        // Total mass.
    double totalDM = 0.0;       // Total diagonal mass.
    double K[3][3];             // Element level stiffness matrix.
    double M[3][3];             // Element level mass matrix.
    double F[3];                // Element level RHS source term vector.
    double x_i[nen], y_i[nen];  // Element-node global coordinates.

    // For the current element:
    //  - zero-initialize M, K, F matrices and vector.
    //  - retrieve the nodal coordinates of the current element.
    for(int i=0; i<nen; i++)
    {
        F[i] = 0.0;
        mesh->getElem(e)->setF(i, 0.0);
        mesh->getElem(e)->setM(i, 0.0);

        for(int j=0; j<nen; j++)
        {
            mesh->getElem(e)->setK(i, j, 0.0);
            K[i][j] = 0.0;
            M[i][j] = 0.0;
        }

        x_i[i] = xyz[mesh->getElem(e)->getConn(i)*nsd + xsd ];
        y_i[i] = xyz[mesh->getElem(e)->getConn(i)*nsd + ysd ];

    }
    
    // Now, calculate the M, K, F matrices.

    // Loop over all Gauss quadrature points.
    for(int p=0; p<nGQP; p++)
    {
        const double zeta = mesh->getME(p)->getPoint(0);        // 1st reference coordinate
        const double eta  = mesh->getME(p)->getPoint(1);        // 2nd reference coordinate
        const double x = calculateRealPosition(zeta,eta,x_i);   // Define x using the isoparametric concept
        const double y = calculateRealPosition(zeta,eta,y_i);   // Define y using the isoparametric concept

        // For the numerical integration, compute a prefactor for the RHS source term vector F. 
        // This prefactor contains:
        // - Gauss integration weights,
        // - the mapping between reference/physical coordinate spaces,
        // - the heat source at the current gauss point.
        double GQweight = mesh->getME(p)->getWeight();
        double detJ = mesh->getElem(e)->getDetJ(p);
        double source = calculateHeatSource(x,y);
        
        const double factor_F = GQweight * detJ * source; // for the heat source contribution.
        
        // Loop over all element nodes (corresponding with weighting functions).
        for(int i=0; i<nen; i++)
        {   
            //Shape function S_i (N_A)
            double S_i = mesh->getME(p)->getS(i);
            double S_i_x = mesh->getElem(e)->getDSdX(p,i);
            double S_i_y = mesh->getElem(e)->getDSdY(p,i);
            // Loop over all element nodes (corresponding with interpolation functions).
            for(int j=0; j<nen; j++)
            {   
                //Shape function S_j (N_B)
                double S_j = mesh->getME(p)->getS(j);
                double S_j_x = mesh->getElem(e)->getDSdX(p,j);
                double S_j_y = mesh->getElem(e)->getDSdY(p,j);
                // Compute element node contributions of the consistent mass matrix.
                M[i][j] += GQweight * S_i * S_j * detJ;

                // Compute element node contributions of the stiffness matrix.
                K[i][j] += D * GQweight * (S_i_x * S_j_x + S_i_y * S_j_y) * detJ;
            }
        
        // Compute element node contributions of the RHS source term vector. Use factor_F and shape function S_i
        F[i] += factor_F * S_i;
        }
    }

    // For the explicit solution, it is necessary to have a diagonal mass matrix, and for this,
    // lumping of the mass matrix is necessary. In order to lump the mass matrix, we first need to
    // calculate the total mass and the total diagonal mass next. Please use the predefined variables:
    // totalM and totalDM
    for(int i=0; i<nen; i++)
    {   
        for(int j=0; j<nen; j++)
        { 
            totalM += M[i][j];
            if(i == j)
                totalDM += M[i][j];
        }
    }

    // Now the diagonal lumping can be done.
    for(int i=0; i<nen; i++)
    {
        for(int j=0; j<nen; j++)
        {   
            if(i == j)
                M[i][j] *= totalM / totalDM; 
            else
                M[i][j] = 0.0;
        }
    }

    //Total mass at each node is accumulated on local node structure:
    for(int i=0; i<nen; i++){
        node = mesh->getElem(e)->getConn(i);
        mesh->getNode(node)->addMass(M[i][i]); //we already have a diagonal matrix. 
    }
    
    // At this point we have the necessary K, M, F matrices as a member of femSolver object.
    // They must be hard copied to the corresponding triElement variables.
    for(int i=0; i<nen; i++)
    {
        mesh->getElem(e)->setF(i, F[i]);
        mesh->getElem(e)->setM(i, M[i][i]);// M[i][i] outside the loop would be sufficient, as M was being diagonalized
        for(int j=0; j<nen; j++)
        {
            mesh->getElem(e)->setK(i, j, K[i][j]);
        }
    }
    // Print the values for element 3921 for debugging purposes.
    if (e==3921)
    {
        for(int i=0; i<nen; i++)
        {
            for(int j=0; j<nen; j++)
            {
                cout << "K: " << K[i][j] << "\t";
            }
            cout << endl;
        }
        for(int i=0; i<nen; i++)
        {
            cout << "M: " << M[i][i] << endl;
        }
        for(int i=0; i<nen; i++)
        {
            cout << "F: " << F[i] << endl;
        }
    }
    return;
}

/***************************************************************************************************
void femSolver::applyDrichletBC()
****************************************************************************************************
if any of the boundary conditions set to Drichlet type:
    visits all partition level nodes:
        determines if any the nodes is on any of the side surfaces.

***************************************************************************************************/
void femSolver::applyDrichletBC()
{
    int const nn = mesh->getNn(); // Number of nodes.
    double * T = mesh->getT(); // Temperature.
    double * xyz = mesh->getXyz(); // Coordinates.
    double x, y, radius; // X coordinate, y coordinate, radius.
    double rOuter = 0.1; // Outer disk radius.
    double temp; // Temporary variable.
    this->nnSolved += 0; // Counter of solved (= non-Dirichlet) nodes.
    //if any of the boundary conditions set to Dirichlet BC:
    if (settings->getBC(1)->getType()==1)
    {   //loop over all nodes.
        for(int i=0; i<nn; i++)
        {   // Calculate the distance of the current node to the center. 
            x = xyz[i*nsd+xsd];
            y = xyz[i*nsd+ysd];
            radius = sqrt(x*x + y*y);
            // Add a threshold for nodes on the outer radius.
            if (abs(radius-rOuter) <= 1E-10)
            { 
                //Set the temperature of the respective node to the value of the Dirichlet BC.
                mesh->getNode(i)->setBCtype(1);
                mesh->getNode(i)->setT(settings->getBC(1)->getValue1()); // 1 equals do outer domain boundary
                // This needs to be synced somehow? Array should contain pointers to individual node temperatures!
                T[i] = settings->getBC(1)->getValue1(); 
            }
            else
            {   
                // Add one to the non-Dirichlet nodes. 
                this->nnSolved++;
            }
        }
    }
    cout << "nnSolved: " << this->nnSolved << endl;

    return;
}

/***************************************************************************************************
* void femSolver::explicitSolver()
***************************************************************************************************
*
**************************************************************************************************/
void femSolver::explicitSolver()
{
    int const nn = mesh->getNn();// Number of nodes.
    int const ne = mesh->getNe();// Number of elements.
    int const nIter = settings->getNIter();//Number of iterations.
    double const dT = settings->getDt();//Time step size.
    double TL[3], MTnewL[3];//Local temperature values of the previous time step and the current time step.
    double * massG = mesh->getMassG(); // Global mass.
    double * MTnew = mesh->getMTnew(); // Temperature at the currently to be evaluated time step. Allocated in triMesh::readMeshFiles()
    double * T = mesh->getT();// Temperature at the previously evaluated time step.
    double massTmp, MTnewTmp;//Temporary variables for mass and M*Tnew.
    double MT; //M*T of the previous time step.
    double Tnew;// Temperature of the currently evaluated time step.
    double partialL2error, globalL2error, initialL2error; // L2 norm of the previous and the current time step, partial (each node), global (all nodes)
    // and initial (time step = 0).
    double* M;          // pointer to element mass matrix.
    double* F;          // pointer to element forcing vector.
    double* K;          // pointer to element stiffness matrix.
    triElement* elem;   // temporary pointer to hold current element.
    triNode* pNode;     // temporary pointer to hold partition nodes.

    for (int iter=0; iter<nIter; iter++)
    {
        // clear RHS MTnew.
        for(int i=0; i<nn; i++){
            MTnew[i] = 0.0; 
        }

        // Evaluate right hand side at element level.
        for(int e=0; e<ne; e++)
        {
            elem = mesh->getElem(e);
            M = elem->getMptr();
            F = elem->getFptr();
            K = elem->getKptr();

            // A bit slower than unrolled alternative
            //prefill TL for subsequent computation
            int Node = 0;
            for(int i=0; i<nen; i++)
            {
                //Get the previous temperature and store it to TL.
                Node = elem->getConn(i);
                TL[i] = T[Node];
            }
            double KT[nen]{};
            for(int i=0; i<nen; i++){
                Node = elem->getConn(i);
                for(int j=0; j<nen; j++)
                    KT[i] += K[i * nen + j] * TL[j];
                
                // Calculate the LHS (M*Tnew) of the matrix equation system (3 entries of the vector per element.)
                MTnewL[i] = M[i] * TL[i] + dT * (F[i] - KT[i]);
                // Save the local values of MT*new (local) (MTNewL) to the global MTNew array.
                MTnew[Node] += MTnewL[i];
            }
            /* for(int i=0; i<nen; i++)
            {
                TL[i] = T[elem->getConn(i)];
            }

            MTnewL[0] = M[0]*TL[0] + dT*(F[0]-(K[0]*TL[0]+K[1]*TL[1]+K[2]*TL[2]));
            MTnewL[1] = M[1]*TL[1] + dT*(F[1]-(K[3]*TL[0]+K[4]*TL[1]+K[5]*TL[2]));
            MTnewL[2] = M[2]*TL[2] + dT*(F[2]-(K[6]*TL[0]+K[7]*TL[1]+K[8]*TL[2]));

            // RHS is accumulated at local nodes
            MTnew[elem->getConn(0)] += MTnewL[0];
            MTnew[elem->getConn(1)] += MTnewL[1];
            MTnew[elem->getConn(2)] += MTnewL[2]; */
        }

        // Evaluate the new temperature on each node on partition level.
        partialL2error = 0.0;
        globalL2error = 0.0;
        // Go over all nodes.
        for(int i=0; i<nn; i++)
        {   //IGNORE Dirichlet nodes.
            pNode = mesh->getNode(i);
            if(pNode->getBCtype() != 1)
            {   // get the mass from the global level. Get M*T from the current timestep and save it to MT. Calculate TNew by inverting the lumped
                // mass matrix.
                MTnewTmp = MTnew[i]; 
                massTmp = massG[i]; //isolate T 1/MT would be inverse of product MT
                
                Tnew = MTnewTmp / massTmp;
                // Calculate the nodal error for the previous and the current timestep.
                partialL2error += pow(Tnew - T[i], 2);
                
                // Set the current temperature to T.
                T[i] = Tnew;
                // Clear MTnew.
                MTnew[i] = 0;
            }
        }
        // Calculate the global error by dividing the partialerror by the number of solved nodes and taking the root.
        globalL2error = sqrt(partialL2error / this->nnSolved);
        
        // Print initial error.
        if(iter==0)
        {
            initialL2error = globalL2error;
            cout << "The initial error is: " << initialL2error << endl;
            cout << "Iter" << '\t' << "Time" << '\t' << "L2 Error" << '\t' << endl;
        }
        //Calculate the relative global error compared to the intial error. 
        globalL2error = globalL2error / initialL2error;
        // Print precision every 1000 time steps. 
        if(iter%1000==0)
        {
            cout << iter << '\t';
            cout << fixed << setprecision(5) << iter*dT << '\t';
            cout << scientific << setprecision(5) << globalL2error << endl;
        }
        // Break the while loop if the global error undergoes a threshold of 1e-7.
        if(globalL2error <= 1.0E-7)
        {
            cout << iter << '\t';
            cout << fixed << setprecision(5) << iter*dT << '\t';
            cout << scientific << setprecision(5) << globalL2error << endl;
            break;
        }
    }
    return;
}


/***************************************************************************************************
* void femSolver::accumulateMass()
***************************************************************************************************
*
**************************************************************************************************/
void femSolver::accumulateMass()
{
    int nn = mesh->getNn();
    double * massG = mesh->getMassG();


    for(int i=0; i<nn; i++)
    {
        massG[i] = mesh->getNode(i)->getMass();
    }

    return;
}

//Position of node?
double femSolver:: calculateRealPosition(const double zeta, const double eta, const double* x) {
    return x[0] + (x[1] - x[0]) * zeta + (x[2] - x[0]) * eta;
}
double femSolver:: calculateHeatSource(const double x, const double y) {
    const double R_12 = settings->getR_1() * settings->getR_1();
    const double r2 = x*x + y*y;
    return (r2 < R_12) ? settings ->getSource() : 0;
}
