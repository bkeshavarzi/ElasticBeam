#include "Q9_Element.h"
using namespace std;
using namespace Eigen;

Q9_Element::Q9_Element()
{
    //ctor
}
Q9_Element::Q9_Element(int i,double thi,vector <Node> obj,ElasticMaterial m)
{
    SetId(i);
    Setth(thi);
    SetNodalObj(obj);
    SetMat(m);
}
int SignFunction(int d)
{
    int res;
    if ((d>0)||(d==0)) res= 1;
    if (d<0) res=-1;
    return res;
}
void Q9_Element::SetId(int i)
{
    id=i;
}
int Q9_Element::GetId(void)
{
    return id;
}
void Q9_Element::Setth(double t)
{
    th=t;
}
double Q9_Element::Getth(void)
{
    return th;
}
void Q9_Element::SetGama(double g)
{
    gama=g;
}
double Q9_Element::GetGama(void)
{
    return gama;
}
void Q9_Element::SetMat(ElasticMaterial m)
{
    mat=m;
    E=mat.GetE();
    v=mat.Getv();
    gama=mat.GetGama();
}
ElasticMaterial Q9_Element::GetMat(void)
{
    return mat;
}
void Q9_Element::SetNodalObj(vector <Node> obj)
{
    NodeObj=obj;
}
vector <Node> Q9_Element::GetNodalObj(void)
{
    return NodeObj;
}
void Q9_Element::SetElemParam(double prop[])
{
    mat.SetE(prop[0]);
    mat.Setv(prop[1]);
    gama=prop[2];
    th=prop[3];
}
double Q9_Element::f_kesi(int kesi_node,double kesi)
{
    double term;
    if (kesi_node==-1) term=kesi*(1-kesi);
    if (kesi_node==0)  term=(1-kesi)*(1+kesi);
    if (kesi_node==1)  term=kesi*(1+kesi);
    return term;
}
double Q9_Element::f_eta(int eta_node,double eta)
{
    double term;
    if (eta_node==-1) term=eta*(1-eta);
    if (eta_node==0)  term=(1-eta)*(1+eta);
    if (eta_node==1)  term=eta*(1+eta);
    return term;
}
double Q9_Element::Calc_ShapeFunction(int kesi_node,int eta_node,double kesi,double eta)
{
    int term1=SignFunction(kesi_node);

    int term2=SignFunction(eta_node);

    double term3=(1/pow(2,abs(kesi_node)+abs(eta_node)));

    return f_kesi(kesi_node,kesi)*f_eta(eta_node,eta)*term1*term2*term3;

}
double Q9_Element::Calc_DiffN(int kesi_node,int eta_node,double kesi,double eta,string str)
{
    double term;

    if ((str=="kesi")&&(kesi_node==-1)) term= (1-2*kesi)*f_eta(eta_node,eta);
    if ((str=="kesi")&&(kesi_node==0))  term= (-1*2*kesi)*f_eta(eta_node,eta);
    if ((str=="kesi")&&(kesi_node==1))  term= (1+2*kesi)*f_eta(eta_node,eta);

    if ((str=="eta")&&(eta_node==-1))   term= (1-2*eta)*f_kesi(kesi_node,kesi);
    if ((str=="eta")&&(eta_node==0))    term= (-1*2*eta)*f_kesi(kesi_node,kesi);
    if ((str=="eta")&&(eta_node==1))    term= (1+2*eta)*f_kesi(kesi_node,kesi);

    double term2=abs(kesi_node)*1.0+abs(eta_node)*1.0;

    return (SignFunction(kesi_node)*SignFunction(eta_node)*(1/pow(2,term2))*term);
}
MatrixXd Q9_Element::Calc_BMatrix(double x_gpt,double y_gpt)
{

    MatrixXd B=MatrixXd::Zero(3,18);
    MatrixXd IJ=MatrixXd::Zero(2,2);
    double kesi,eta;

    for (int inode=0;inode<9;inode++)
    {
        kesi=(LCord(0,inode)-xc)/a;
        eta=(LCord(1,inode)-yc)/b;
        IJ=CalcInvJacobian(x_gpt,y_gpt);
        B(0,2*inode)=IJ(0,0)*Calc_DiffN(kesi,eta,x_gpt,y_gpt,"kesi")+IJ(0,1)*Calc_DiffN(kesi,eta,x_gpt,y_gpt,"eta");
        B(1,2*inode+1)=IJ(1,0)*Calc_DiffN(kesi,eta,x_gpt,y_gpt,"kesi")+IJ(1,1)*Calc_DiffN(kesi,eta,x_gpt,y_gpt,"eta");
        B(2,2*inode)=B(1,2*inode+1);
        B(2,2*inode+1)=B(0,2*inode);
    }
    return B;
}
void Q9_Element::SetDMatrix(string str)
{
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
        D(0,0)=E/(1-pow(v,2.0));
        D(0,1)=E*v/(1-pow(v,2.0));
        D(1,0)=D(0,1);
        D(1,1)=D(0,0);
        D(2,2)=E/(1+v);
    }
}
void Q9_Element::Setlocalcord(void)
{
    vector <double> Cord;

    for (int inode=0;inode<9;inode++)
    {
        Cord=NodeObj[inode].GetCord();
        LCord(0,inode)=Cord[0];
        LCord(1,inode)=Cord[1];
    }
    double max_x=LCord.block(0,0,1,9).maxCoeff();
    double min_x=LCord.block(0,0,1,9).minCoeff();
    double max_y=LCord.block(1,0,1,9).maxCoeff();
    double min_y=LCord.block(1,0,1,9).minCoeff();
    a=0.5*(max_x-min_x);
    b=0.5*(max_y-min_y);
    xc=0.5*(max_x+min_x);
    yc=0.5*(max_y+min_y);
}
MatrixXd Q9_Element::CalcJacobian(double kesi,double eta)
{
    double kesinode;
    double etanode;
    MatrixXd J=MatrixXd::Zero(2,2);

    for (int inode=0;inode<9;inode++)
    {
        kesinode=(LCord(0,inode)-xc)/a;
        etanode=(LCord(1,inode)-yc)/b;
        J(0,0)+=LCord(0,inode)*Calc_DiffN(kesinode,etanode,kesi,eta,"kesi");
        J(0,1)+=LCord(1,inode)*Calc_DiffN(kesinode,etanode,kesi,eta,"kesi");
        J(1,0)+=LCord(0,inode)*Calc_DiffN(kesinode,etanode,kesi,eta,"eta");
        J(1,1)+=LCord(1,inode)*Calc_DiffN(kesinode,etanode,kesi,eta,"eta");
    }
    return J;
}
MatrixXd Q9_Element::CalcInvJacobian(double kesi,double eta)
{
    return (CalcJacobian(kesi,eta).inverse());
}
double Q9_Element::CalcDetJacobian(double kesi,double eta)
{
    return (CalcJacobian(kesi,eta).determinant());
}
MatrixXd Q9_Element::GetLocalCord()
{
    return LCord;
}
MatrixXd Q9_Element::Calc_LSM()
{
    MatrixXd B=MatrixXd::Zero(3,18);
    double detJ;

    for (int igpt=0;igpt<3;igpt++)
    {
        for (int jpt=0;jpt<3;jpt++)
        {
            B=Calc_BMatrix(gpt[igpt],gpt[jpt]);
            detJ=CalcDetJacobian(gpt[igpt],gpt[jpt]);
            LSM+=wgpt[igpt]*(B.transpose())*D*B*detJ;
        }
    }

    LSM=((LSM.array())*th).matrix();

}
void Q9_Element::SetU(MatrixXd Ue)
{
    U=Ue;
}
MatrixXd Q9_Element::GetU()
{
    return U;
}
void Q9_Element::Set_Ep()
{
    Ep.block(0,0,3,1)=Calc_BMatrix(gpt[0],gpt[0])*U;
    Ep.block(0,1,3,1)=Calc_BMatrix(gpt[1],gpt[0])*U;
    Ep.block(0,2,3,1)=Calc_BMatrix(gpt[2],gpt[0])*U;

    Ep.block(0,3,3,1)=Calc_BMatrix(gpt[0],gpt[1])*U;
    Ep.block(0,4,3,1)=Calc_BMatrix(gpt[1],gpt[1])*U;
    Ep.block(0,5,3,1)=Calc_BMatrix(gpt[2],gpt[1])*U;

    Ep.block(0,6,3,1)=Calc_BMatrix(gpt[0],gpt[2])*U;
    Ep.block(0,7,3,1)=Calc_BMatrix(gpt[1],gpt[2])*U;
    Ep.block(0,8,3,1)=Calc_BMatrix(gpt[2],gpt[2])*U;
}
MatrixXd Q9_Element::Get_Ep()
{
    return Ep;
}
void Q9_Element::Set_Sigma()
{
    Sigma.block(0,0,3,1)=D*Ep.block(0,0,3,1);
    Sigma.block(0,1,3,1)=D*Ep.block(0,1,3,1);
    Sigma.block(0,2,3,1)=D*Ep.block(0,2,3,1);

    Sigma.block(0,3,3,1)=D*Ep.block(0,3,3,1);
    Sigma.block(0,4,3,1)=D*Ep.block(0,4,3,1);
    Sigma.block(0,5,3,1)=D*Ep.block(0,5,3,1);

    Sigma.block(0,6,3,1)=D*Ep.block(0,6,3,1);
    Sigma.block(0,7,3,1)=D*Ep.block(0,7,3,1);
    Sigma.block(0,8,3,1)=D*Ep.block(0,8,3,1);

}
MatrixXd Q9_Element::Get_Sigma()
{
    return Sigma;
}
void Q9_Element::Set_PSigma(string str)
{
    MatrixXd Sigma_Ave=0.5*(Sigma.block(0,0,2,9).colwise().sum()); //n x 1
    MatrixXd R=pow(pow(0.5*(Sigma.block(0,0,1,9).array()-Sigma.block(1,0,1,9).array()),2.0)+pow(Sigma.block(2,0,1,9).array(),2.0),0.5).matrix();
    PSigma.block(0,0,1,9)=Sigma_Ave+R;
    PSigma.block(1,0,1,9)=Sigma_Ave-R;
    if (str=="pstrain") PSigma.block(2,0,1,9)=(v*(PSigma.block(0,0,1,9)+PSigma.block(1,0,1,9)).array()).matrix();
    if (str=="pstress") PSigma.block(2,0,1,9)=MatrixXd::Zero(1,9);
}
MatrixXd Q9_Element::Get_PSigma()
{
    return PSigma;
}
void Q9_Element::Set_PStrain(string str)
{
    MatrixXd Strain_Ave=0.5*(Ep.block(0,0,2,9).colwise().sum()); //n x 1
    MatrixXd R=pow(pow(0.5*(Ep.block(0,0,1,9).array()-Ep.block(1,0,1,9).array()),2.0)+pow(Ep.block(2,0,1,9).array(),2.0),0.5).matrix();
    PStrain.block(0,0,1,9)=Strain_Ave+R;
    PStrain.block(1,0,1,9)=Strain_Ave-R;
    if (str=="pstress") PStrain.block(2,0,1,9)=((v/E)*(PSigma.block(0,0,1,9)+PSigma.block(1,0,1,9)).array()).matrix();
    if (str=="pstrain") PStrain.block(2,0,1,9)=MatrixXd::Zero(1,9);
}
MatrixXd Q9_Element::Get_PStrain()
{
    return PStrain;
}
Q9_Element::~Q9_Element()
{
    //dtor
}
