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
void Q4_Element::CalcJacobian(double kesi,double eta)
{
    J(0,0)=localcord(0,0)*Calc_DiffN(-1,-1,kesi,eta,"kesi")+localcord(0,1)*Calc_DiffN(1,-1,kesi,eta,"kesi")+localcord(0,2)*Calc_DiffN(1,1,kesi,eta,"kesi")+localcord(0,3)*Calc_DiffN(-1,1,kesi,eta,"kesi");
    J(0,1)=localcord(1,0)*Calc_DiffN(-1,-1,kesi,eta,"kesi")+localcord(1,1)*Calc_DiffN(1,-1,kesi,eta,"kesi")+localcord(1,2)*Calc_DiffN(1,1,kesi,eta,"kesi")+localcord(1,3)*Calc_DiffN(-1,1,kesi,eta,"kesi");
    J(1,0)=localcord(0,0)*Calc_DiffN(-1,-1,kesi,eta,"eta")+localcord(0,1)*Calc_DiffN(1,-1,kesi,eta,"eta")+localcord(0,2)*Calc_DiffN(1,1,kesi,eta,"eta")+localcord(0,3)*Calc_DiffN(-1,1,kesi,eta,"eta");
    J(1,1)=localcord(1,0)*Calc_DiffN(-1,-1,kesi,eta,"eta")+localcord(1,1)*Calc_DiffN(1,-1,kesi,eta,"eta")+localcord(1,2)*Calc_DiffN(1,1,kesi,eta,"eta")+localcord(1,3)*Calc_DiffN(-1,1,kesi,eta,"eta");
}
void Q4_Element::CalcInvJacobian(double kesi,double eta)
{
    CalcJacobian(double kesi,double eta);
    IJ=J.inverse();
}
void Q4_Element::CalcDetJacobian(double kesi,double eta)
{
    CalcJacobian(double kesi,double eta);
    detJ=J.determinant();
}
double Q4_Element::Calc_DiffN(int kesi_node,int eta_node,double kesi,double eta,string str)
{
    if (str=="kesi") return (0.25*kesi_node*(1+eta_node*eta));
    if (str=="eta")  return (0.25*eta_node*(1+kesi_node*kesi));
}
MatrixXd Q4_Element::Calc_BMatrix(double x_gpt,double y_gpt)
{
    CalcInvJacobian(x_gpt,y_gpt);
    B(0,0)=IJ(0,0)*Calc_DiffN(-1,-1,x_gpt,y_gpt,"kesi")+IJ(0,1)*Calc_DiffN(-1,-1,x_gpt,y_gpt,"eta");
    B(1,1)=IJ(1,0)*Calc_DiffN(-1,-1,x_gpt,y_gpt,"kesi")+IJ(1,1)*Calc_DiffN(-1,-1,x_gpt,y_gpt,"eta");
    B(2,0)=B(1,1);B(2,1)=B(0,0);
    B(0,2)=IJ(0,0)*Calc_DiffN(1,-1,x_gpt,y_gpt,"kesi")+IJ(0,1)*Calc_DiffN(1,-1,x_gpt,y_gpt,"eta");
    B(1,3)=IJ(1,0)*Calc_DiffN(1,-1,x_gpt,y_gpt,"kesi")+IJ(1,1)*Calc_DiffN(1,-1,x_gpt,y_gpt,"eta");
    B(2,2)=B(1,3);B(2,3)=B(0,2);
    B(0,4)=IJ(0,0)*Calc_DiffN(1,1,x_gpt,y_gpt,"kesi")+IJ(0,1)*Calc_DiffN(1,1,x_gpt,y_gpt,"eta");
    B(1,5)=IJ(1,0)*Calc_DiffN(1,1,x_gpt,y_gpt,"kesi")+IJ(1,1)*Calc_DiffN(1,1,x_gpt,y_gpt,"eta");
    B(2,4)=B(1,5);B(2,5)=B(0,4);
    B(0,6)=IJ(0,0)*Calc_DiffN(-1,1,x_gpt,y_gpt,"kesi")+IJ(0,1)*Calc_DiffN(-1,1,x_gpt,y_gpt,"eta");
    B(1,7)=IJ(1,0)*Calc_DiffN(-1,1,x_gpt,y_gpt,"kesi")+IJ(1,1)*Calc_DiffN(-1,1,x_gpt,y_gpt,"eta");
    B(2,6)=B(1,7);B(2,7)=B(0,6);
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
    double xmin,xmax,ymin,ymax,xsum,ysum;
    cord=NodeObj[0].GetCord();
    xmin=cord[0];xmax=cord[0];ymin=cord[1];ymax=cord[1];xsum=cord[0];ysum=cord[1];
    localcord(0,0)=cord[0];localcord(1,0)=cord[1];

    for (int inode=1;inode<4;inode++)
    {
        cord=NodeObj[inode].GetCord();
        xmin=(xmin>cord[0]?cord[0]:xmin);
        xmax=(xmax>cord[0]?xmax:cord[0]);
        ymin=(ymin>cord[1]?cord[1]:ymin);
        ymax=(ymax>cord[1]?ymax:cord[1]);
        xsum+=cord[0];
        ysum+=cord[1];
    }
    a=0.5*(xmax-xmin);
    b=0.5*(ymax-ymin);
    xsum=xsum/4;
    ysum=ysum/4;
    localcord.block(0,0,1,4)=(((localcord.block(0,0,1,4)).array()-xsum)/a).matrix();
    localcord.block(1,0,1,4)=(((localcord.block(1,0,1,4)).array()-ysum)/b).matrix();
}
MatrixXd Q4_Element::Getlocalcord()
{
    return localcord;
}
void Q4_Element::Calc_LSM()
{
    MatrixXd TempB=MatrixXd::Zero(3,8);
    TempB=Calc_BMatrix(gpt[0],gpt[0]);
    CalcDetJacobian(gpt[0],gpt[0]);
    LSM=gpt[0]*(TempB.transpose())*D*TempB*detJ;

    TempB=Calc_BMatrix(gpt[1],gpt[0]);
    CalcDetJacobian(gpt[1],gpt[0]);
    LSM+=gpt[1]*(TempB.transpose())*D*TempB*detJ;

    TempB=Calc_BMatrix(gpt[1],gpt[1]);
    CalcDetJacobian(gpt[1],gpt[1]);
    LSM+=gpt[1]*(TempB.transpose())*D*TempB*detJ;

    TempB=Calc_BMatrix(gpt[0],gpt[1]);
    CalcDetJacobian(gpt[0],gpt[1]);
    LSM+=gpt[0]*(TempB.transpose())*D*TempB;

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
