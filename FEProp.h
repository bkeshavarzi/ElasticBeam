#ifndef FEPROP_H_INCLUDED
#define FEPROP_H_INCLUDED
#include <iostream>

#include <Eigen/Dense>
#include "Node.h"
#include "CSE_Element.h"
#include "Q4_Element.h"

using namespace std;
using namespace Eigen;

MatrixXd AssembleForceVector(vector <Node>);
MatrixXd AssembleStiffnessMatrix_CSE(vector <Node> NV,vector <CSE_Element>);
MatrixXd AssembleStiffnessMatrix_Q4(vector <Q4_Element>);

#endif // FEPROP_H_INCLUDED
