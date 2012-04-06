#include "CSE_Element.h"

CSE_Element::CSE_Element()
{
    //vector <Node> ObjNode=CSE_Element_v[0].GetNodeObj();
    //CSE_Element obj(1,0.02,ObjNode,mat);
    //obj.SetDMatrix_PlaneStress();
    //obj.SetLSM();
}
CSE_Element::CSE_Element(int i,double thi,vector <Node> obj,ElasticMaterial m)
{
    SetId(i);
    Setth(thi);
    SetNodeObj(obj[0],obj[1],obj[2]);
    SetMat(m);
    SetCMatrix();
    SetA();
    SetBMatrix();
}
void CSE_Element::SetId(int i)
{
    id=i;
}
int CSE_Element::GetId(void)
{
    return id;
}
void CSE_Element::Setth(double thi)
{
    th=thi;
}
double CSE_Element::Getth(void)
{
    return th;
}
void CSE_Element::SetGama(double g)
{
    gama=g;
}
double CSE_Element::GetGama(void)
{
    return gama;
}
void CSE_Element::SetMat(ElasticMaterial m)
{
    mat=m;
    E=m.GetE();
    v=m.Getv();
    gama=mat.GetGama();
}
ElasticMaterial CSE_Element::GetMat(void)
{
    return mat;
}
void CSE_Element::SetNodeObj(Node obj1, Node obj2, Node obj3)
{
    NodeObj.push_back(obj1);
    NodeObj.push_back(obj2);
    NodeObj.push_back(obj3);
}
vector <Node> CSE_Element::GetNodeObj(void)
{
    return NodeObj;
}
void CSE_Element::SetDof()
{
    *ldof=2*(NodeObj[0].GetId()-1);
    *(ldof+1)=2*(NodeObj[0].GetId()-1)+1;
    *(ldof+2)=2*(NodeObj[1].GetId()-1);
    *(ldof+3)=2*(NodeObj[1].GetId()-1)+1;
    *(ldof+4)=2*(NodeObj[2].GetId()-1);
    *(ldof+5)=2*(NodeObj[2].GetId()-1)+1;
}
void CSE_Element::SetCMatrix()
{
    vector <double> Temp_Cord;
    Temp_Cord=NodeObj[0].GetCord();
    C(0,0)=1;C(0,1)=Temp_Cord[0];C(0,2)=Temp_Cord[1];
    Temp_Cord=NodeObj[1].GetCord();
    C(1,0)=1;C(1,1)=Temp_Cord[0];C(1,2)=Temp_Cord[1];
    Temp_Cord=NodeObj[2].GetCord();
    C(2,0)=1;C(2,1)=Temp_Cord[0];C(2,2)=Temp_Cord[1];
}
MatrixXd CSE_Element::GetCMatrix()
{
    return C;
}
void CSE_Element::SetA()
{
    A=0.5*(C.determinant());
}
double CSE_Element::GetA()
{
    return A;
}
void CSE_Element::SetBMatrix()
{
    MatrixXd Cinv=C.inverse();
    B(0,0)=Cinv(1,0);B(1,1)=Cinv(2,0);B(2,0)=Cinv(2,0);B(2,1)=Cinv(1,0);
    B(0,2)=Cinv(1,1);B(1,3)=Cinv(2,1);B(2,2)=Cinv(2,1);B(2,3)=Cinv(1,1);
    B(0,4)=Cinv(1,2);B(1,5)=Cinv(2,2);B(2,4)=Cinv(2,2);B(2,5)=Cinv(1,2);
}
MatrixXd CSE_Element::GetBMatrix(void)
{
    return B;
}
void CSE_Element::SetDMatrix_PlaneStress()
{
    D(0,0)=(E/(1-pow(v,2.0)));
    D(0,1)=E*v/(1-pow(v,2.0));
    D(1,0)=D(0,1);
    D(1,1)=D(0,0);
    D(2,2)=E/(1+v);
}
void CSE_Element::SetDMatrix_PlaneStrain()
{
    double factor=E/((1+v)*(1-2*v));
    D(0,0)=factor*(1-v);
    D(0,1)=factor*v;
    D(1,0)=D(0,1);
    D(1,1)=D(0,0);
    D(2,2)=factor*(1-2*v);
}
MatrixXd CSE_Element::GetDMatrix()
{
    return D;
}
void CSE_Element::SetLSM()
{
    LSM=B.transpose()*D*B*th*A;
}
MatrixXd CSE_Element::GetLSM()
{
    return LSM;
}
void CSE_Element::SetU(MatrixXd Ue)
{
    U=Ue;
}
void CSE_Element::CalculateStrainTensor()
{
    Ep=B*U;
}
MatrixXd CSE_Element::GetStrainTensor()
{
    return Ep;
}
void CSE_Element::CalculateStressTensor()
{
    Sigma=D*Ep;
}
MatrixXd CSE_Element::GetStressTensor()
{
    return Sigma;
}
void CSE_Element::CalculatePrinStress(string ptype)
{
    double sigma_ave,R,sigmaZ;

    if (ptype=="plane strain") double sigmaZ=v*(Sigma(0)+Sigma(1));
    if (ptype=="plane stress") double sigmaZ=0;

    sigma_ave=0.5*(Sigma(0)+Sigma(1));
    R=pow(pow(0.5*(Sigma(0)-Sigma(1)),2)+pow(Sigma(2),2),0.5);

    PSigma(0)=sigma_ave-R;
    PSigma(1)=sigma_ave+R;
    PSigma(2)=sigmaZ;
}
MatrixXd CSE_Element::GetPrinStress()
{
    return PSigma;
}
void CSE_Element::CalculatePrinStrain(string ptype)
{
    double strain_ave,R,strainZ;

    if (ptype=="plane strain") strainZ=0;
    if (ptype=="plane stress") strainZ=(-v/E)*(Sigma(0)+Sigma(1));

    strain_ave=0.5*(Ep(0)+Ep(1));
    R=pow(pow(0.5*(Ep(0)-Ep(1)),2)+pow(Ep(2),2),0.5);

    PStrain(0)=strain_ave+R;
    PStrain(1)=strain_ave-R;
    PStrain(2)=strainZ;

}
MatrixXd CSE_Element::GetPrinStrain()
{
    return PStrain;
}
CSE_Element::~CSE_Element()
{
    //dtor
}
