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

    cout << "Total number of nodes are: " <<NV.size()<<endl;
    cout << "Total number of elements are: " <<EV.size()<<endl;

    MatrixXd FDOF=BoundryCondition(NV);
    MatrixXd FV=AssembleForceVector(NV);

    MatrixXd KG=AssembleStiffnessMatrix_CSE(NV,EV);
    MatrixXd FKG=CondenseStiffnessMatrix(NV,KG,FDOF);
    MatrixXd FFV=CondenseForceVector(NV,FV,FDOF);

    cout << "Number of Restrained DOF are : " << 2*NV.size()-FDOF.sum()<<endl;
    //cout << "Number of elements in force vector are : " << FV.minCoeff()<<endl;
    //cout << "Number of rows and columns in stiffness are : " << KG.rows()<< "\t" << KG.cols() << endl;
    cout << "Number of free dof in stiffness are : " << FKG.rows()<< "\t" << FKG.cols() << endl;
    //cout << "Number of free dof in force are : " << FFV.rows()<<endl;

    //MatrixXd FU=Solve_CSE(NV,EV,FFV,FKG,FDOF);
    cout << (FKG.squaredNorm())*((FKG.inverse()).squaredNorm()) <<endl;
    //cout << FU.minCoeff() <<endl;
    //cout << FU.maxCoeff() <<endl;

    //WriteOutPutFile_CSE(FU,NV,EV);
}
