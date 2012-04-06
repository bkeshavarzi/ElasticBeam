#include "Q4_Element.h"
using namespace std;
using namespace Eigen;

Q4_Element::Q4_Element()
{
    //ctor
}
Q4_Element::Q4_Element(int i,double thi,double g, vector <Node> obj)
{
    SetId(i);
    Setth(thi);
    SetGama(g);
    SetNodalObj(obj);
}
void Q4_Element::SetId(int i)
{
    id=i;
}
int Q4_Element::GetId(void)
{
    return id;
}
void Q4_Element::Setth(double t)
{
    th=t;
}
double Q4_Element::Getth(void)
{
    return th;
}
void Q4_Element::SetGama(double g)
{
    gama=g;
}
double Q4_Element::GetGama(void)
{
    return gama;
}
void Q4_Element::SetMat(ElasticMaterial m)
{
    mat=m;
    E=mat.GetE();
    v=mat.Getv();
    gama=mat.GetGama();
}
ElasticMaterial Q4_Element::GetMat(void)
{
    return mat;
}
void Q4_Element::SetNodalObj(vector <Node> obj)
{
    NodeObj=obj;
}
double Q4_Element::Calc_ShapeFunction(int kesi_node,int eta_node,double kesi,double eta)
{
    return 0.25*(1+kesi_node*kesi)*(1+eta_node*eta);
}
double Q4_Element::Calc_DiffN(int kesi_node,int eta_node,double kesi,double eta,string str)
{
    if (str=="kesi") return (0.25*kesi_node*(1+eta_node*eta));
    if (str=="eta")  return (0.25*eta_node*(1+kesi_node*kesi));
}
MatrixXd Q4_Element::Calc_BMatrix(double x_gpt,double y_gpt)
{
    B(0,0)=Calc_DiffN(-1,-1,x_gpt,y_gpt,"kesi");B(1,1)=Calc_DiffN(-1,-1,x_gpt,y_gpt,"eta");B(2,0)=B(1,1);B(2,1)=B(0,0);
    B(0,2)=Calc_DiffN(1,-1,x_gpt,y_gpt,"kesi");B(1,3)=Calc_DiffN(1,-1,x_gpt,y_gpt,"eta");B(2,2)=B(1,3);B(2,3)=B(0,2);
    B(0,4)=Calc_DiffN(1,1,x_gpt,y_gpt,"kesi");B(1,5)=Calc_DiffN(1,1,x_gpt,y_gpt,"eta");B(2,4)=B(1,5);B(2,5)=B(0,4);
    B(0,6)=Calc_DiffN(-1,1,x_gpt,y_gpt,"kesi");B(1,7)=Calc_DiffN(-1,1,x_gpt,y_gpt,"eta");B(2,6)=B(1,5);B(2,7)=B(0,6);

    return B;
}
void Q4_Element::SetDMatrix(string str)
{
    MatrixXd D=MatrixXd::Zero(3,3);
    if (str=="pstrain")
    {
        double factor=E/((1+v)*(1-2*v));
        D(0,0)=factor*(1-v);
        D(0,1)=factor*v;
        D(1,0)=D(0,1);
        D(1,1)=D(0,0);
        D(2,2)=factor*(1-2*v);
    }
    if (str=="pstress")
    {
        D(0,0)=(E/(1-pow(v,2)));
        D(0,1)=E*v/(1-pow(v,2));
        D(1,0)=D(0,1);
        D(1,1)=D(0,0);
        D(2,2)=E/(1+v);
    }
}
void Q4_Element::Calc_LSM()
{
    MatrixXd TempB=MatrixXd::Zero(3,8);
    TempB=Calc_BMatrix(gpt[0],gpt[0]);
    LSM=TempB.transpose()*D*TempB;
    TempB=Calc_BMatrix(gpt[1],gpt[0]);
    LSM+=TempB.transpose()*D*TempB;
    TempB=Calc_BMatrix(gpt[1],gpt[1]);
    LSM+=TempB.transpose()*D*TempB;
    TempB=Calc_BMatrix(gpt[0],gpt[1]);
    LSM+=TempB.transpose()*D*TempB;
}
MatrixXd Q4_Element::Get_LSM()
{
    return LSM;
}
void Q4_Element::SetU(MatrixXd Ue)
{
    U=Ue;
}
MatrixXd Q4_Element::Get_Ep()
{
    Ep.block(0,0,3,1)=Calc_BMatrix(gpt[0],gpt[0])*U;
    Ep.block(0,1,3,1)=Calc_BMatrix(gpt[1],gpt[0])*U;
    Ep.block(0,2,3,1)=Calc_BMatrix(gpt[1],gpt[1])*U;
    Ep.block(0,3,3,1)=Calc_BMatrix(gpt[0],gpt[1])*U;

    return Ep;
}
MatrixXd Q4_Element::Get_Sigma()
{
    Sigma.block(0,0,3,1)=D*Ep.block(0,0,3,1);
    Sigma.block(0,1,3,1)=D*Ep.block(0,1,3,1);
    Sigma.block(0,2,3,1)=D*Ep.block(0,2,3,1);
    Sigma.block(0,3,3,1)=D*Ep.block(0,3,3,1);

    return Sigma;
}
MatrixXd Q4_Element::Get_PSigma(string str)
{
    MatrixXd Sigma_Ave=0.5*(Sigma.block(0,0,2,4).colwise().sum()); //n x 1
    MatrixXd R=pow(pow(0.5*(Sigma.block(0,0,1,4).array()-Sigma.block(1,0,1,4).array()),2.0)+pow(Sigma.block(2,0,1,4).array(),2.0),0.5).matrix();
    PSigma.block(0,0,1,4)=Sigma_Ave+R;
    PSigma.block(1,0,1,4)=Sigma_Ave-R;
    if (str=="pstrain") PSigma.block(2,0,1,4)=(v*(PSigma.block(0,0,1,4)+PSigma.block(1,0,1,4)).array()).matrix();
    if (str=="pstress") PSigma.block(2,0,1,4)=MatrixXd::Zero(1,4);

    return PSigma;
}
MatrixXd Q4_Element::Get_PStrain(string str)
{
    MatrixXd Strain_Ave=0.5*(Ep.block(0,0,2,4).colwise().sum()); //n x 1
    MatrixXd R=pow(pow(0.5*(Ep.block(0,0,1,4).array()-Ep.block(1,0,1,4).array()),2.0)+pow(Ep.block(2,0,1,4).array(),2.0),0.5).matrix();
    PStrain.block(0,0,1,4)=Strain_Ave+R;
    PStrain.block(1,0,1,4)=Strain_Ave-R;
    if (str=="pstress") PStrain.block(2,0,1,4)=((v/E)*(PSigma.block(0,0,1,4)+PSigma.block(1,0,1,4)).array()).matrix();
    if (str=="pstrain") PStrain.block(2,0,1,4)=MatrixXd::Zero(1,4);

    return PStrain;
}
Q4_Element::~Q4_Element()
{
    //dtor
}
