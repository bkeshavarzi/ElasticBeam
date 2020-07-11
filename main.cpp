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


//CSE element works very well, Stiffness matrix checked.


int main()
{
    ElasticMaterial mat(2e9,0.3,2500);
    ElasticMaterial & m=mat;
    vector <Node> NV=ReadNodeFile("Q8_Node.txt");
    vector <Q8_Element> EV=ReadQ8Element("Q8_Element.txt",m,0.02,NV);
    //cout << NV.size() <<"\t" << EV.size() <<endl;
    //cout << m.GetE()<< "\t"<< m.Getv() <<"\t"<< m.GetGama() <<endl;
    //cout << EV[0].GetId() << "\t" << EV[0].Getth() << "\t" << EV[0].GetGama() << "\t" << endl;
    //vector <Node> obj=EV[0].GetNodalObj();
    //cout << obj[0].GetId() << "\t" << obj[1].GetId() <<"\t"<< obj[2].GetId() <<"\t"<<obj[3].GetId()<<"\t"<<obj[4].GetId()<<"\t"<<obj[5].GetId()<<"\t"<<obj[6].GetId()<<"\t"<<obj[7].GetId()<<endl;
     //cout << EV[0].GetLocalCord() << endl;
    //cout << EV[0].Calc_ShapeFunction(-1,-1,5,5)<<endl;
    //cout << EV[0].Calc_DiffN(-1,-1,5,5,"kesi")<<endl;
    //EV[0].SetA();
    //cout << EV[0].GetA() <<endl;
    //cout << EV[0].Calc_DiffN(0,1,-.5,-0.5,"kesi") <<endl;
     //cout << EV[0].Calc_BMatrix(-1/sqrt(3),-1/sqrt(3)) <<endl;
     //cout << EV[0].Calc_DiffN(1,1,-0.5,-0.5,"kesi") <<endl;
     cout << EV[0].CalcJacobian(-1/sqrt(3),1/sqrt(3)) <<endl;
    //cout << EV[0].CalcInvJacobian(5,5) <<endl;
    //cout << EV[0].CalcDetJacobian(5,5);
    //cout << EV[0].Calc_LSM()<<endl;
    //cout << EV[0].Calc_LSM() <<endl;
}
