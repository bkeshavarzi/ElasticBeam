#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "Node.h"
#include "CSE_Element.h"
#include "Q4_Element.h"
#include "Q8_Element.h"
#include "Q9_Element.h"
#include "TextFiles.h"

using namespace std;
using namespace Eigen;

//Make node class...Done
//Mesh structures...DONE
//Make element class(Q4...DONE, Q8...DONE, Q9...DONE, CSE...DONE)...DONE
//Read node file...DONE
//Read element file...DONE
//A code in Python to show mesh and results in contour mode, try to link it with C++
//Make stiffness matrix
//Make force vector
//Make free dof vector
//Think about node load
//Calculate displacement vector
//Calculate stress for each element
//Calculate prin. stress for each stress

int main()
{
    //ElasticMaterial mat(2e9,0.3,0);
    //vector <Node> NV=ReadNodeFile("Q9_Node.txt");
    //vector <Q9_Element> EV=ReadQ9Element("Q9_Element.txt",mat,0.02,NV);
    //vector <Node> nv=EV[0].GetNodalObj();
    //EV[0].Setlocalcord();
    //cout << EV[0].Calc_DiffN(-1,-1,-1*sqrt(0.6),-1*sqrt(0.6),"eta") << endl;
    //MatrixXd B =EV[0].Calc_BMatrix(-1*sqrt(0.6),-1*sqrt(0.6));
    //cout << B.block(0,0,3,2) <<endl;
}
