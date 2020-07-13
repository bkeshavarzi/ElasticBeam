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
//Make stiffness matrix...DONE
//Make force vector...DONE
//Make free dof vector...DONE
//Think about node load...DONE
//Calculate displacement vector...DONE
//Calculate stress for each element...DONE
//Calculate prin. stress for each stress...DONE


int main()
{
    ElasticMaterial mat(2e9,0.3,2500);
    ElasticMaterial & m=mat;
}
