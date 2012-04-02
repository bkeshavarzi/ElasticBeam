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
    B(0,0)=Calc_DiffN(-1,-1,x_gpt,y_gpt,"kesi");B(1,1)=Calc_DiffN(-1,-1,x_gpt,y_gpt,"eta");B(2,0)=B(1,1);B(2,1)=B(0,0);
    B(0,2)=Calc_DiffN(0,-1,x_gpt,y_gpt,"kesi");B(1,3)=Calc_DiffN(0,-1,x_gpt,y_gpt,"eta");B(2,2)=B(1,3);B(2,3)=B(0,2);
    B(0,4)=Calc_DiffN(1,-1,x_gpt,y_gpt,"kesi");B(1,5)=Calc_DiffN(1,-1,x_gpt,y_gpt,"eta");B(2,4)=B(1,5);B(2,5)=B(0,4);
    B(0,6)=Calc_DiffN(1,0,x_gpt,y_gpt,"kesi");B(1,7)=Calc_DiffN(1,0,x_gpt,y_gpt,"eta");B(2,6)=B(1,7);B(2,7)=B(0,6);
    B(0,8)=Calc_DiffN(0,0,x_gpt,y_gpt,"kesi");B(1,9)=Calc_DiffN(0,0,x_gpt,y_gpt,"eta");B(2,8)=B(1,9);B(2,9)=B(0,8);
    B(0,10)=Calc_DiffN(-1,0,x_gpt,y_gpt,"kesi");B(1,11)=Calc_DiffN(-1,0,x_gpt,y_gpt,"eta");B(2,10)=B(1,11);B(2,11)=B(0,10);
    B(0,12)=Calc_DiffN(-1,1,x_gpt,y_gpt,"kesi");B(1,13)=Calc_DiffN(-1,1,x_gpt,y_gpt,"eta");B(2,12)=B(1,13);B(2,13)=B(0,12);
    B(0,14)=Calc_DiffN(0,1,x_gpt,y_gpt,"kesi");B(1,15)=Calc_DiffN(0,1,x_gpt,y_gpt,"eta");B(2,14)=B(1,15);B(2,15)=B(0,14);
    B(0,16)=Calc_DiffN(1,1,x_gpt,y_gpt,"kesi");B(1,17)=Calc_DiffN(1,1,x_gpt,y_gpt,"eta");B(2,16)=B(1,17);B(2,17)=B(0,16);

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
    vector <double> cord;
    double xmin,xmax,ymin,ymax,xsum,ysum;
    cord=NodeObj[0].GetCord();
    xmin=cord[0];xmax=cord[0];ymin=cord[1];ymax=cord[1];xsum=cord[0];ysum=cord[1];
    localcord(0,0)=cord[0];localcord(1,0)=cord[1];

    for (int inode=1;inode<9;inode++)
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
    localcord.block(0,0,1,9)=(((localcord.block(0,0,1,9)).array()-xsum)/a).matrix();
    localcord.block(1,0,1,9)=(((localcord.block(1,0,1,9)).array()-ysum)/b).matrix();
}
MatrixXd Q9_Element::Getlocalcord()
{
    return localcord;
}
void Q9_Element::Calc_LSM()
{
    LSM=wgpt[0]*(Calc_BMatrix(gpt[0],gpt[0]).transpose()*D*Calc_BMatrix(gpt[0],gpt[0])+Calc_BMatrix(gpt[1],gpt[0]).transpose()*D*Calc_BMatrix(gpt[1],gpt[0])+Calc_BMatrix(gpt[2],gpt[0]).transpose()*D*Calc_BMatrix(gpt[2],gpt[0]));
    LSM=LSM+wgpt[1]*(Calc_BMatrix(gpt[0],gpt[1]).transpose()*D*Calc_BMatrix(gpt[0],gpt[1])+Calc_BMatrix(gpt[1],gpt[1]).transpose()*D*Calc_BMatrix(gpt[1],gpt[1])+Calc_BMatrix(gpt[2],gpt[1]).transpose()*D*Calc_BMatrix(gpt[2],gpt[1]));
    LSM=LSM+wgpt[2]*(Calc_BMatrix(gpt[0],gpt[2]).transpose()*D*Calc_BMatrix(gpt[0],gpt[2])+Calc_BMatrix(gpt[1],gpt[2]).transpose()*D*Calc_BMatrix(gpt[1],gpt[2])+Calc_BMatrix(gpt[2],gpt[2]).transpose()*D*Calc_BMatrix(gpt[2],gpt[2]));
    LSM=th*LSM;
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
