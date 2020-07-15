#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "Node.h"
#include "CSE_Element.h"
#include "Q4_Element.h"
#include "Q8_Element.h"
#include "Q9_Element.h"
#include "TextFiles.h"
#include "FEProp.h"

using namespace std;
using namespace Eigen;


int main()
{
    ElasticMaterial mat(2e9,0.3,2500);
    ElasticMaterial & m=mat;

    vector <Node> NV=ReadNodeFile("CSE_Node.txt");
    vector <CSE_Element> EV=ReadCSEElement("CSE_Element.txt",m,0.002,NV);

    //cout << "Total number of nodes are: " <<NV.size()<<endl;
    //cout << "Total number of elements are: " <<EV.size()<<endl;


    MatrixXd FDOF=BoundryCondition(NV);
    MatrixXd FV=AssembleForceVector(NV);
    MatrixXd KG=AssembleStiffnessMatrix_CSE(NV,EV);
    MatrixXd FKG=CondenseStiffnessMatrix_CSE(NV,KG,FDOF);
    MatrixXd FFV=CondenseForceVector_CSE(NV,FV,FDOF);

    MatrixXd UG=MatrixXd::Zero(FDOF.sum(),1);
    MatrixXd UT=MatrixXd::Zero(2*NV.size(),1);
    MatrixXd U=MatrixXd::Zero(6,1);

    Solve_CSE(NV,EV,FFV,FKG,FDOF);
    //WriteOutPutFile_CSE(U,NV,EV);

}
