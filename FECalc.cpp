#include "FEProp.h"

MatrixXd AssembleForceVector(vector <Node> NV)
{
    MatrixXd FV=MatrixXd::Zero(2*NV.size(),1);
    vector <double> temp;
    temp.push_back(0.0);temp.push_back(0.0);
    for (int inode=0;inode<NV.size();inode++)
    {
        temp=NV[inode].GetForce();
        FV(2*inode)=temp[0];
        FV(2*inode+1)=temp[1];
    }
    return FV;
}

MatrixXd AssembleStiffnessMatrix_CSE(vector <Node> NV,vector <CSE_Element> EV)
{
    MatrixXd KG=MatrixXd::Zero(2*NV.size(),2*NV.size());
    MatrixXd KE=MatrixXd::Zero(6,6);
    vector <Node> NodeObj;
    MatrixXd LDOF=MatrixXd::Zero(6,1);

    for (int ielem=0; ielem<EV.size(); ielem++)
    {
        KE=EV[ielem].Calc_LSM();
        NodeObj=EV[ielem].GetNodeObj();

        for (int inode=0;inode<3;inode++)
        {
            LDOF(2*inode)=2*(NodeObj[inode].GetId()-1);
            LDOF(2*inode+1)=2*(NodeObj[inode].GetId()-1)+1;
        }
        for (int idof=0;idof<6;idof++)
        {
            for (int jdof=0;jdof<6;jdof++)
            {
                KG(LDOF(idof),LDOF(jdof))+=KE(idof,jdof);
            }
        }
    }
    return KG;
}

MatrixXd AssembleStiffnessMatrix_Q4(vector <Node> NV,vector <Q4_Element> EV)
{
    MatrixXd KG=MatrixXd::Zero(2*NV.size(),2*NV.size());
    MatrixXd KE=MatrixXd::Zero(8,8);
    vector <Node> NodeObj;
    MatrixXd LDOF=MatrixXd::Zero(8,1);

    for (int ielem=0; ielem<EV.size(); ielem++)
    {

        NodeObj=EV[ielem].GetNodalObj();
        KE=EV[ielem].Calc_LSM();
        for (int inode=0;inode<4;inode++)
        {
            LDOF(2*inode)=2*(NodeObj[inode].GetId()-1);
            LDOF(2*inode+1)=2*(NodeObj[inode].GetId()-1)+1;
        }

        for (int idof=0;idof<8;idof++)
        {
            for (int jdof=0;jdof<8;jdof++)
            {
                KG(LDOF(idof),LDOF(jdof))+=KE(idof,jdof);
            }
        }
    }
    return KG;
}

MatrixXd AssembleStiffnessMatrix_Q8(vector <Node> NV,vector <Q8_Element> EV)
{
    MatrixXd KG=MatrixXd::Zero(2*NV.size(),2*NV.size());
    MatrixXd KE=MatrixXd::Zero(16,16);
    vector <Node> NodeObj;
    MatrixXd LDOF=MatrixXd::Zero(16,1);

    for (int ielem=0; ielem<EV.size(); ielem++)
    {
        NodeObj=EV[ielem].GetNodalObj();
        KE=EV[ielem].Calc_LSM();
        for (int inode=0;inode<8;inode++)
        {
            LDOF(2*inode)=2*(NodeObj[inode].GetId()-1);
            LDOF(2*inode+1)=2*(NodeObj[inode].GetId()-1)+1;
        }
        for (int idof=0;idof<16;idof++)
        {
            for (int jdof=0;jdof<16;jdof++)
            {
                KG(LDOF(idof),LDOF(jdof))+=KE(idof,jdof);
            }
        }
    }

}

MatrixXd AssembleStiffnessMatrix_Q9(vector <Node> NV,vector <Q9_Element> EV)
{
    MatrixXd KG=MatrixXd::Zero(2*NV.size(),2*NV.size());
    MatrixXd KE=MatrixXd::Zero(18,18);
    vector <Node> NodeObj;
    MatrixXd LDOF=MatrixXd::Zero(18,1);

    for (int ielem=0; ielem<EV.size(); ielem++)
    {
        NodeObj=EV[ielem].GetNodalObj();
        KE=EV[ielem].Calc_LSM();
        for (int inode=0;inode<9;inode++)
        {
            LDOF(2*inode)=2*(NodeObj[inode].GetId()-1);
            LDOF(2*inode+1)=2*(NodeObj[inode].GetId()-1)+1;
        }
        for (int idof=0;idof<18;idof++)
        {
            for (int jdof=0;jdof<18;jdof++)
            {
                KG(LDOF(idof),LDOF(jdof))+=KE(idof,jdof);
            }
        }
    }

}
