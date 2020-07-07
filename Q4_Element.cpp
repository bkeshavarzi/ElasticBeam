#include "Q4_Element.h"
using namespace std;
using namespace Eigen;

Q4_Element::Q4_Element()
{
    //ctor
}
Q4_Element::Q4_Element(int i,double thi, vector <Node> obj,ElasticMaterial mat)
{
    SetId(i);
    Setth(thi);
    SetMat(mat);
    SetNodalObj(obj);
    Setlocalcord();
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
vector <Node> Q4_Element::GetNodalObj(void)
{
    return NodeObj;
}
double Q4_Element::Calc_ShapeFunction(int kesi_node,int eta_node,double kesi,double eta)
{
    return 0.25*(1+kesi_node*kesi)*(1+eta_node*eta);
}
MatrixXd Q4_Element::CalcJacobian(double kesi,double eta)
{
    MatrixXd J=MatrixXd::Zero(2,2);
    int kesi_node,eta_node;

    for (int inode=0;inode<4;inode++)
    {
        kesi_node=(LCord(0,inode)-xc)/a;
        eta_node=(LCord(1,inode)-yc)/b;
        J(0,0)+=LCord(0,inode)*Calc_DiffN(kesi_node,eta_node,kesi,eta,"kesi");
        J(0,1)+=LCord(1,inode)*Calc_DiffN(kesi_node,eta_node,kesi,eta,"kesi");
        J(1,0)+=LCord(0,inode)*Calc_DiffN(kesi_node,eta_node,kesi,eta,"eta");
        J(1,1)+=LCord(1,inode)*Calc_DiffN(kesi_node,eta_node,kesi,eta,"eta");
    }
    return J;
}
MatrixXd Q4_Element::CalcInvJacobian(double kesi,double eta)
{
    MatrixXd I,IJ=MatrixXd::Zero(2,2);
    J=CalcJacobian(double kesi,double eta);
    IJ=J.inverse();
    return IJ;
}
double Q4_Element::CalcDetJacobian(double kesi,double eta)
{
    MatrixXd J=MatrixXd::Zero(2,2);
    J=CalcJacobian(double kesi,double eta);
    double detJ=J.determinant();
    return detJ;
}
double Q4_Element::Calc_DiffN(int kesi_node,int eta_node,double kesi,double eta,string str)
{
    if (str=="kesi") return (0.25*kesi_node*(1+eta_node*eta));
    if (str=="eta")  return (0.25*eta_node*(1+kesi_node*kesi));
}
MatrixXd Q4_Element::Calc_BMatrix(double x_gpt,double y_gpt)
{
    MatrixXd IJ=CalcInvJacobian(x_gpt,y_gpt);
    int kesi_node,eta_node;
    MatrixXd B=MatrixXd::Zero(3,8);

    for (int inode=0;inode<4;inode++)
    {
        kesi_node=(LCord(0,ionde)-xc)/a;
        eta_node=(LCord(1,ionde)-yc)/b;
        B(0,2*inode)=IJ(0,0)*Calc_DiffN(kesi_node,eta_node,x_gpt,y_gpt,"kesi")+IJ(0,1)*Calc_DiffN(kesi_node,eta_node,x_gpt,y_gpt,"eta");
        B(1,2*inode+1)=IJ(1,0)*Calc_DiffN(kesi_node,eta_node,x_gpt,y_gpt,"kesi")+IJ(1,1)*Calc_DiffN(kesi_node,eta_node,x_gpt,y_gpt,"eta");
        B(2,2*inode)=B(1,2*inode+1);
        B(2,2*inode+1)=B(0,2*inode);
    }
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
void Q4_Element::Setlocalcord(void)
{
    vector <double> cord;

    for (int inode=0;inode<4;inode++)
    {
        cord=NodeObj[inode].GetCord();
        LCord(0,inode)=cord[0];
        LCord(1,inode)=cord[1];
    }
    a=0.5*(LCord.block(0,0,1,4).maxCoeff()-LCord.block(0,0,1,4).minCoeff());
    b=0.5*(LCord.block(1,0,1,4).maxCoeff()-LCord.block(1,0,1,4).minCoeff());
    xc=(LCord.block(0,0,1,4).sum())/4;
    yc=(LCord.block(0,1,1,4).sum())/4;
}
void Q4_Element::Calc_LSM()
{
    double detJ;
    MatrixXd B=MatrixXd::Zero(3,8);
    for (int igpt=0;igpt <2;igpt++)
    {
        for (int jgpt=0;jgpt<2;jgpt++)
        {
            B=Calc_BMatrix(gpt[igpt],gpt[jgpt]);
            detJ=CalcDetJacobian(gpt[igpt],gpt[jgpt]);
            LSM+=detJ*(B.transpose())*D*B;
        }
    }
    LSM=((LSM.array())*th).matrix();
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
    MatrixXd Sigma_Ave=0.5*(Sigma.block(0,0,2,4).colwise().sum()); //1 x n
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
