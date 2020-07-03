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

    MatrixXd B=MatrixXd::Zero(3,16);
    MatrixXd IJ=MatrixXd::Zeros(2,2);
    double kesi,eta;

    for (int inode=0;inode<8;inode++)
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
void Q8_Element::SetLocalCord(void)
{
    vector <double> Cord;

    for (int inode=0;inode<8;inode++)
    {
        Cord=NodeObj[inode].GetCord();
        LCord(0,inode)=Cord[0];
        LCord(1,inode)=Cord[1];
    }
    double max_x=LCord.block(0,0,1,8).maxCoeff();
    double min_x=LCord.block(0,0,1,8).minCoeff();
    double max_y=LCord.block(1,0,1,8).maxCoeff();
    double min_y=LCord.block(1,0,1,8).minCoeff();
    a=0.5*(max_x-min_x);
    b=0.5*(max_y-min_y);
    xc=0.5*(max_x+min_x);
    yc=0.5*(max_y+min_y);
}
MatrixXd Q8_Element::CalcJacobian(double kesi,double eta)
{
    double kesinode;
    double etanode;
    MatrixXd J=MatrixXd::Zero(2,2);

    for (int inode=0;inode<8;inode++)
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
void Q8_Element::CalcInvJacobian(double kesi,double eta)
{
    return (CalcJacobian(kesi,eta).inverse());
}
double Q8_Element::CalcDetJacobian(double kesi,double eta)
{
    return (CalcJacobian(kesi,eta).determinant());
}
void Q8_Element::Calc_LSM()
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
