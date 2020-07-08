#ifndef FEPROP_H_INCLUDED
#define FEPROP_H_INCLUDED
#include <iostream>

#include <Eigen/Dense>
#include "Node.h"
#include "CSE_Element.h"
#include "Q4_Element.h"
#include "Q8_Element.h"
#include "Q9_Element.h"

using namespace std;
using namespace Eigen;

MatrixXd AssembleForceVector(vector <Node>);
MatrixXd AssembleStiffnessMatrix_CSE(vector <Node> NV,vector <CSE_Element>);
MatrixXd AssembleStiffnessMatrix_Q4(vector <Node>,vector <Q4_Element>);
MatrixXd AssembleStiffnessMatrix_Q8(vector <Node>,vector <Q8_Element>);
MatrixXd AssembleStiffnessMatrix_Q9(vector <Node>,vector <Q9_Element>);

#endif // FEPROP_H_INCLUDED
