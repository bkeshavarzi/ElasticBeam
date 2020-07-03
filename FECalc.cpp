#include "FEProp.h"

MatrixXd AssembleForceVector(vector <Node> NV)
{
    MatrixXd FV=MatrixXd::Zero(2*NV.size(),1);
    vector <double> temp;
    temp.push_back(0.0);temp.push_back(0.0);
    for (int inode=0;inode<NV.size();inode++)
    {
        temp.NV[inode].GetForce();
        FV(2*inode)=temp[0];
        FV(2*inode+1)=temp[1];
    }
    return FV;
}

MatrixXd AssembleStiffnessMatrix_CSE(vector <Node> NV,vector <CSE_Element> EV)
{
    MatrixXd KG=MatrixXd::Zero(2*NV.size(),2*NV.size());
    MatrixXd KE=MatrixXd::Zero(6,6);
    int * LDOF;
    for (int ielem=0; ielem<EV.size(); ielem++)
    {
        EV[ielem].SetFECalc("planestress");
        KE=EV[ielem].SetLSM();
        LDOF=EV[ielem].ldof;
        for (int idof=0;idof<6;idof++)
        {
            for (int jdof=0;jdof<6;jdof++)
            {
                KG(LDOF[idof],LDOF[jdof])+=KE(idof,jdof);
            }
        }
    }
    return KG;
}

MatrixXd AssembleStiffnessMatrix_Q4(vector <Node> NV,vector <Q4_Element> EV)
{
    MatrixXd KG=MatrixXd::Zero(2*NV.size(),2*NV.size());
    MatrixXd KE=MatrixXd::Zero(8,8);
    int LDOF[8]=[0,0,0,0,0,0,0,0];
    vector <Node> EN;

    for (int ielem=0; ielem<EV.size(); ielem++)
    {
        EV[ielem].SetDMatrix("pstress");
        EV[ielem].Calc_LSM();
        KE=EV[ielem].Get_LSM();
        EN=EV[ielem].GetNodalObj();

        LDOF[0]=2*(EN[0].GetId()-1);
        LDOF[1]=2*(EN[0].GetId()-1)+1;
        LDOF[2]=2*(EN[1].GetId()-1);
        LDOF[3]=2*(EN[1].GetId()-1)+1;
        LDOF[4]=2*(EN[2].GetId()-1);
        LDOF[5]=2*(EN[2].GetId()-1)+1;
        LDOF[6]=2*(EN[3].GetId()-1);
        LDOF[7]=2*(EN[3].GetId()-1)+1;

        for (int idof=0;idof<8;idof++)
        {
            for (int jdof=0;jdof<8;jdof++)
            {
                KG(LDOF[idof],LDOF[jdof])+=KE(idof,jdof);
            }
        }
    }
    return KG;
}
