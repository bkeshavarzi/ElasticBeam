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
int SignFunction(double d)
{
    int res;
    if ((d<0)||(d==0)) res= 1;
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
double Q9_Element::Calc_ShapeFunction(int kesi_node,int eta_node,double kesi,double eta)
{
    double term1,term2;

    if (kesi_node==-1) term1=kesi*(kesi-1);
    if (kesi_node==0)  term1=(kesi-1)*(kesi+1);
    if (kesi_node==1)  term1=kesi*(kesi+1);

    if (eta_node==-1)  term2=eta*(eta-1);
    if (eta_node==0)   term2=(eta+1)*(eta-1);
    if (eta_node==1)   term2=(eta+1)*(eta);

    cout << term1 << "\t" << term2 << "\t" <<SignFunction(kesi_node)*SignFunction(eta_node) << "\t" <<(1/pow(2,abs(kesi_node)+abs(eta_node)))<<endl;

    return SignFunction(kesi_node)*SignFunction(eta_node)*(1/pow(2,abs(kesi_node)+abs(eta_node)))*term1*term2;
}
double Q9_Element::Calc_DiffN(int kesi_node,int eta_node,double kesi,double eta,string str)
{
    double term;
    if ((str=="kesi")&&(kesi_node==-1)) term= (2*kesi-1)*(eta_node==-1?eta*(eta-1):1)*(eta_node==0?(eta+1)*(eta-1):1)*(eta_node==1?(eta+1)*(eta):1);
    if ((str=="kesi")&&(kesi_node==0))  term= (2*kesi)*(eta_node==-1?eta*(eta-1):1)*(eta_node==0?(eta+1)*(eta-1):1)*(eta_node==1?(eta+1)*(eta):1);
    if ((str=="kesi")&&(kesi_node==1))  term= (2*kesi+1)*(eta_node==-1?eta*(eta-1):1)*(eta_node==0?(eta+1)*(eta-1):1)*(eta_node==1?(eta+1)*(eta):1);
    if ((str=="eta")&&(eta_node==-1))   term= (2*eta-1)*(kesi_node==-1?kesi*(kesi-1):1)*(kesi_node==0?(kesi+1)*(kesi-1):1)*(kesi_node==1?(kesi+1)*(kesi):1);
    if ((str=="eta")&&(eta_node==0))    term= (2*eta)*(kesi_node==-1?kesi*(kesi-1):1)*(kesi_node==0?(kesi+1)*(kesi-1):1)*(kesi_node==1?(kesi+1)*(kesi):1);
    if ((str=="eta")&&(eta_node==1))    term= (2*eta+1)*(kesi_node==-1?kesi*(kesi-1):1)*(kesi_node==0?(kesi+1)*(kesi-1):1)*(kesi_node==1?(kesi+1)*(kesi):1);

    return SignFunction(kesi_node)*SignFunction(eta_node)*(1/pow(2,abs(kesi_node)+abs(eta_node)))*term;
}
MatrixXd Q9_Element::Calc_BMatrix(double x_gpt,double y_gpt)
{

    MatrixXd B=MatrixXd::Zero(3,18);
    MatrixXd IJ=MatrixXd::Zeros(2,2);
    double kesi,eta;

    for (int inode=0;inode<9;inode++)
    {
        kesi=(LCord(0,inode)-xc)/a;
        eta=(LCord(1,inode)-yc)/b);
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
double Q9_Element::CalcDetJacobian(double,double)
{
    return (CalcJacobian(kesi,eta).determinant());
}
MatrixXd Q9_Element::Getlocalcord()
{
    return localcord;
}
void Q9_Element::Calc_LSM()
{
    for (int igpt=0;igpt<3;igpt++)
    {
        for (int jpt=0;jpt<3;jpt++)
        {
            LSM+=wgpt[igpt]*(Calc_BMatrix(gpt[igpt],gpt[jpt]).transpose())*D*Calc_BMatrix(gpt[igpt],gpt[jpt])*CalcDetJacobian(gpt[igpt],gpt[jpt]);
        }
    }

    LSM=LSM*th;

}
MatrixXd Q9_Element::Get_LSM()
{
    return LSM;
}
void Q9_Element::SetU(MatrixXd Ue)
{
    U=Ue;
}
MatrixXd Q9_Element::Get_Ep()
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

    return Ep;
}
MatrixXd Q9_Element::Get_Sigma()
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

    return Sigma;
}
MatrixXd Q9_Element::Get_PSigma(string str)
{
    MatrixXd Sigma_Ave=0.5*(Sigma.block(0,0,2,9).colwise().sum()); //n x 1
    MatrixXd R=pow(pow(0.5*(Sigma.block(0,0,1,9).array()-Sigma.block(1,0,1,9).array()),2.0)+pow(Sigma.block(2,0,1,9).array(),2.0),0.5).matrix();
    PSigma.block(0,0,1,9)=Sigma_Ave+R;
    PSigma.block(1,0,1,9)=Sigma_Ave-R;
    if (str=="pstrain") PSigma.block(2,0,1,9)=(v*(PSigma.block(0,0,1,9)+PSigma.block(1,0,1,9)).array()).matrix();
    if (str=="pstress") PSigma.block(2,0,1,9)=MatrixXd::Zero(1,9);

    return PSigma;
}
MatrixXd Q9_Element::Get_PStrain(string str)
{
    MatrixXd Strain_Ave=0.5*(Ep.block(0,0,2,9).colwise().sum()); //n x 1
    MatrixXd R=pow(pow(0.5*(Ep.block(0,0,1,9).array()-Ep.block(1,0,1,9).array()),2.0)+pow(Ep.block(2,0,1,9).array(),2.0),0.5).matrix();
    PStrain.block(0,0,1,9)=Strain_Ave+R;
    PStrain.block(1,0,1,9)=Strain_Ave-R;
    if (str=="pstress") PStrain.block(2,0,1,9)=((v/E)*(PSigma.block(0,0,1,9)+PSigma.block(1,0,1,9)).array()).matrix();
    if (str=="pstrain") PStrain.block(2,0,1,9)=MatrixXd::Zero(1,9);

    return PStrain;
}
Q9_Element::~Q9_Element()
{
    //dtor
}
