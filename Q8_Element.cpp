#include "Q8_Element.h"

Q8_Element::Q8_Element()
{
    //ctor
}
Q8_Element::Q8_Element(int i,double thi,vector <Node> obj)
{
    SetId(i);
    Setth(thi);
    SetNodalObj(obj);
}
void Q8_Element::SetId(int i)
{
    id=i;
}
int Q8_Element::GetId(void)
{
    return id;
}
void Q8_Element::Setth(double t)
{
    th=t;
}
double Q8_Element::Getth(void)
{
    return th;
}
void Q8_Element::SetGama(double g)
{
    gama=g;
}
double Q8_Element::GetGama(void)
{
    return gama;
}
void Q8_Element::SetMat(ElasticMaterial m)
{
    mat=m;
    E=mat.GetE();
    v=mat.Getv();
    gama=mat.GetGama();
}
ElasticMaterial Q8_Element::GetMat(void)
{
    return mat;
}
void Q8_Element::SetNodalObj(vector <Node> obj)
{
    NodeObj=obj;
}
vector <Node> Q8_Element:: GetNodalObj(void)
{
    return NodeObj;
}
double Q8_Element::Calc_ShapeFunction(int kesi_node,int eta_node,double kesi,double eta) //double Calc_ShapeFunction(int,int,double,double);
{
    double term=0;

    if ((kesi_node==-1)&&(eta_node==-1))   term=(Q4_Element::Calc_ShapeFunction(-1,-1,kesi,eta))-0.5*(Q9_Element::Calc_ShapeFunction(0,-1,kesi,eta))-0.5*(Q9_Element::Calc_ShapeFunction(-1,0,kesi,eta));
    if ((kesi_node==0)&&(eta_node==-1))    term=0.5*(kesi+1)*(kesi-1)*(eta-1);
    if ((kesi_node==1)&&(eta_node==-1))    term=Q4_Element::Calc_ShapeFunction(1,-1,kesi,eta)-0.5*Q9_Element::Calc_ShapeFunction(0,-1,kesi,eta)-0.5*Q9_Element::Calc_ShapeFunction(1,0,kesi,eta);
    if ((kesi_node==1)&&(eta_node==0))     term=-0.5*(eta-1)*(kesi+1)*(eta+1);
    if ((kesi_node==-1)&&(eta_node==0))    term=0.5*(eta-1)*(kesi-1)*(eta+1);
    if ((kesi_node==-1)&&(eta_node==1))    term=Q4_Element::Calc_ShapeFunction(-1,1,kesi,eta)-0.5*Q9_Element::Calc_ShapeFunction(-1,0,kesi,eta)-0.5*Q9_Element::Calc_ShapeFunction(0,1,kesi,eta);
    if ((kesi_node==0)&&(eta_node==1))     term=-0.5*(kesi+1)*(kesi-1)*(eta+1);
    if ((kesi_node==1)&&(eta_node==1))     term=Q4_Element::Calc_ShapeFunction(1,1,kesi,eta)-0.5*Q9_Element::Calc_ShapeFunction(0,1,kesi,eta)-0.5*Q9_Element::Calc_ShapeFunction(1,0,kesi,eta);

    return term;
}
double Q8_Element::Calc_DiffN(int kesi_node,int eta_node,double kesi,double eta,string str)
{
    double term;

    if ((str=="kesi")&&(kesi_node==-1)&&(eta_node==-1)) term=Q4_Element::Calc_DiffN(-1,-1,kesi,eta,"kesi")-0.5*Q9_Element::Calc_DiffN(0,-1,kesi,eta,"kesi")-0.5*Q9_Element::Calc_DiffN(-1,0,kesi,eta,"kesi");
    if ((str=="eta")&&(kesi_node==-1)&&(eta_node==-1))  term=Q4_Element::Calc_DiffN(-1,-1,kesi,eta,"eta")-0.5*Q9_Element::Calc_DiffN(0,-1,kesi,eta,"eta")-0.5*Q9_Element::Calc_DiffN(-1,0,kesi,eta,"eta");
    if ((str=="kesi")&&(kesi_node==0)&&(eta_node==-1))  term=kesi*(eta-1);
    if ((str=="eta")&&(kesi_node==0)&&(eta_node==-1))   term=0.5*(pow(kesi,2)-1);
    if ((str=="kesi")&&(kesi_node==1)&&(eta_node==-1))  term=Q4_Element::Calc_DiffN(1,-1,kesi,eta,"kesi")-0.5*Q9_Element::Calc_DiffN(0,-1,kesi,eta,"kesi")-0.5*Q9_Element::Calc_DiffN(1,0,kesi,eta,"kesi");
    if ((str=="eta")&&(kesi_node==1)&&(eta_node==-1))   term=Q4_Element::Calc_DiffN(1,-1,kesi,eta,"eta")-0.5*Q9_Element::Calc_DiffN(0,-1,kesi,eta,"eta")-0.5*Q9_Element::Calc_DiffN(1,0,kesi,eta,"eta");
    if ((str=="kesi")&&(kesi_node==1)&&(eta_node==0))   term=-0.5*(pow(eta,2)-1);
    if ((str=="eta")&&(kesi_node==1)&&(eta_node==0))    term=-1*eta*(kesi+1);
    if ((str=="kesi")&&(kesi_node==-1)&&(eta_node==0))  term=0.5*(pow(eta,2)-1);
    if ((str=="eta")&&(kesi_node==-1)&&(eta_node==0))   term=eta*(kesi-1);
    if ((str=="kesi")&&(kesi_node==-1)&&(eta_node==1))  term=Q4_Element::Calc_DiffN(-1,1,kesi,eta,"kesi")-0.5*Q9_Element::Calc_DiffN(-1,0,kesi,eta,"kesi")-0.5*Q9_Element::Calc_DiffN(0,1,kesi,eta,"kesi");
    if ((str=="eta")&&(kesi_node==-1)&&(eta_node==1))   term=Q4_Element::Calc_DiffN(-1,1,kesi,eta,"eta")-0.5*Q9_Element::Calc_DiffN(-1,0,kesi,eta,"eta")-0.5*Q9_Element::Calc_DiffN(0,1,kesi,eta,"eta");
    if ((str=="kesi")&&(kesi_node==0)&&(eta_node==1))   term=-0.5*kesi*(eta+1);
    if ((str=="eta")&&(kesi_node==0)&&(eta_node==1))    term=-0.5*(pow(kesi,2)-1);
    if ((str=="kesi")&&(kesi_node==1)&&(eta_node==1))   term=Q4_Element::Calc_DiffN(1,1,kesi,eta,"kesi")-0.5*Q9_Element::Calc_DiffN(0,1,kesi,eta,"kesi")-0.5*Q9_Element::Calc_DiffN(1,0,kesi,eta,"kesi");
    if ((str=="eta")&&(kesi_node==1)&&(eta_node==1))    term=Q4_Element::Calc_DiffN(1,1,kesi,eta,"eta")-0.5*Q9_Element::Calc_DiffN(0,1,kesi,eta,"eta")-0.5*Q9_Element::Calc_DiffN(1,0,kesi,eta,"eta");

    return term;
}
MatrixXd Q8_Element::Calc_BMatrix(double x_gpt,double y_gpt)
{
    B(0,0)=Calc_DiffN(-1,-1,x_gpt,y_gpt,"kesi");B(1,1)=Calc_DiffN(-1,-1,x_gpt,y_gpt,"eta");B(2,0)=B(1,1);B(2,1)=B(0,0);
    B(0,2)=Calc_DiffN(0,-1,x_gpt,y_gpt,"kesi");B(1,3)=Calc_DiffN(0,-1,x_gpt,y_gpt,"eta");B(2,2)=B(1,3);B(2,3)=B(0,2);
    B(0,4)=Calc_DiffN(1,-1,x_gpt,y_gpt,"kesi");B(1,5)=Calc_DiffN(1,-1,x_gpt,y_gpt,"eta");B(2,4)=B(1,5);B(2,5)=B(0,4);
    B(0,6)=Calc_DiffN(1,0,x_gpt,y_gpt,"kesi");B(1,7)=Calc_DiffN(1,0,x_gpt,y_gpt,"eta");B(2,6)=B(1,7);B(2,7)=B(0,6);
    B(0,8)=Calc_DiffN(0,0,x_gpt,y_gpt,"kesi");B(1,9)=Calc_DiffN(0,0,x_gpt,y_gpt,"eta");B(2,6)=B(1,7);B(2,7)=B(0,6);
    B(0,10)=Calc_DiffN(-1,0,x_gpt,y_gpt,"kesi");B(1,11)=Calc_DiffN(-1,0,x_gpt,y_gpt,"eta");B(2,8)=B(1,9);B(2,9)=B(0,8);
    B(0,12)=Calc_DiffN(-1,1,x_gpt,y_gpt,"kesi");B(1,11)=Calc_DiffN(-1,1,x_gpt,y_gpt,"eta");B(2,10)=B(1,11);B(2,11)=B(0,10);
    B(0,14)=Calc_DiffN(0,1,x_gpt,y_gpt,"kesi");B(1,13)=Calc_DiffN(0,1,x_gpt,y_gpt,"eta");B(2,12)=B(1,13);B(2,13)=B(0,12);
    B(0,16)=Calc_DiffN(1,1,x_gpt,y_gpt,"kesi");B(1,15)=Calc_DiffN(1,1,x_gpt,y_gpt,"eta");B(2,14)=B(1,15);B(2,15)=B(0,14);

    return B;
}
void Q8_Element::SetDMatrix(string str)
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
void Q8_Element::Calc_LSM()
{
    LSM=wgpt[0]*(Calc_BMatrix(gpt[0],gpt[0]).transpose()*D*Calc_BMatrix(gpt[0],gpt[0])+Calc_BMatrix(gpt[1],gpt[0]).transpose()*D*Calc_BMatrix(gpt[1],gpt[0])+Calc_BMatrix(gpt[2],gpt[0]).transpose()*D*Calc_BMatrix(gpt[2],gpt[0]));
    LSM=LSM+wgpt[1]*(Calc_BMatrix(gpt[0],gpt[1]).transpose()*D*Calc_BMatrix(gpt[0],gpt[1])+Calc_BMatrix(gpt[1],gpt[1]).transpose()*D*Calc_BMatrix(gpt[1],gpt[1])+Calc_BMatrix(gpt[2],gpt[1]).transpose()*D*Calc_BMatrix(gpt[2],gpt[1]));
    LSM=LSM+wgpt[2]*(Calc_BMatrix(gpt[0],gpt[2]).transpose()*D*Calc_BMatrix(gpt[0],gpt[2])+Calc_BMatrix(gpt[1],gpt[2]).transpose()*D*Calc_BMatrix(gpt[1],gpt[2])+Calc_BMatrix(gpt[2],gpt[2]).transpose()*D*Calc_BMatrix(gpt[2],gpt[2]));
    LSM=th*LSM;
}
MatrixXd Q8_Element::Get_LSM()
{
    return LSM;
}
void Q8_Element::SetU(MatrixXd Ue)
{
    U=Ue;
}
MatrixXd Q8_Element::Get_Ep()
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
MatrixXd Q8_Element::Get_Sigma()
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
MatrixXd Q8_Element::Get_PSigma(string str)
{
    MatrixXd Sigma_Ave=0.5*(Sigma.block(0,0,2,9).colwise().sum()); //n x 1
    MatrixXd R=pow(pow(0.5*(Sigma.block(0,0,1,9).array()-Sigma.block(1,0,1,9).array()),2.0)+pow(Sigma.block(2,0,1,9).array(),2.0),0.5).matrix();
    PSigma.block(0,0,1,9)=Sigma_Ave+R;
    PSigma.block(1,0,1,9)=Sigma_Ave-R;
    if (str=="pstrain") PSigma.block(2,0,1,9)=(v*(PSigma.block(0,0,1,9)+PSigma.block(1,0,1,9)).array()).matrix();
    if (str=="pstress") PSigma.block(2,0,1,9)=MatrixXd::Zero(1,9);

    return PSigma;
}
MatrixXd Q8_Element::Get_PStrain(string str)
{
    MatrixXd Strain_Ave=0.5*(Ep.block(0,0,2,9).colwise().sum()); //n x 1
    MatrixXd R=pow(pow(0.5*(Ep.block(0,0,1,9).array()-Ep.block(1,0,1,9).array()),2.0)+pow(Ep.block(2,0,1,9).array(),2.0),0.5).matrix();
    PStrain.block(0,0,1,9)=Strain_Ave+R;
    PStrain.block(1,0,1,9)=Strain_Ave-R;
    if (str=="pstress") PStrain.block(2,0,1,9)=((v/E)*(PSigma.block(0,0,1,9)+PSigma.block(1,0,1,9)).array()).matrix();
    if (str=="pstrain") PStrain.block(2,0,1,9)=MatrixXd::Zero(1,9);

    return PStrain;
}
Q8_Element::~Q8_Element()
{
    //dtor
}
