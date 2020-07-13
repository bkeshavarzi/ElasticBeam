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

MatrixXd BoundryCondition(vector <Node>);
MatrixXd AssembleForceVector(vector <Node>);
MatrixXd AssembleStiffnessMatrix_CSE(vector <Node> NV,vector <CSE_Element>);
MatrixXd AssembleStiffnessMatrix_Q4(vector <Node>,vector <Q4_Element>);
MatrixXd AssembleStiffnessMatrix_Q8(vector <Node>,vector <Q8_Element>);
MatrixXd AssembleStiffnessMatrix_Q9(vector <Node>,vector <Q9_Element>);
MatrixXd CondenseStiffnessMatrix(vector <Node>,MatrixXd,MatrixXd);
MatrixXd Solve_CSE(vector <Node>,vector <CSE_Element>,MatrixXd,MatrixXd,MatrixXd);
MatrixXd Solve_Q4(vector <Node>,vector <Q4_Element>,MatrixXd,MatrixXd,MatrixXd);
MatrixXd Solve_Q9(vector <Node>,vector <Q4_Element>,MatrixXd,MatrixXd,MatrixXd);
#endif // FEPROP_H_INCLUDED
