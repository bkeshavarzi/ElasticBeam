#include "FEProp.h"


MatrixXd BoundryCondition(vector <Node> NV)
{
    vector <double> temp;
    MatrixXd FDOF=MatrixXd::Zero(2*NV.size(),1);

    for (int inode=0;inode<NV.size();inode++)
    {
        temp=NV[inode].GetCord();
        // 0=restraiend, 1=free;
        if (temp[0]!=0)
        {
            FDOF(2*inode)=1;
            FDOF(2*inode+1)=1;
        }
    }
    return FDOF;
}
MatrixXd AssembleForceVector(vector <Node> NV)
{
    MatrixXd FV=MatrixXd::Zero(2*NV.size(),1);
    vector <double> temp;

    for (int inode=0;inode<NV.size();inode++)
    {
        temp=NV[inode].GetCord();
        if ((temp[0]==3)&&(temp[0]==0.5))
        {
            FV(2*inode)=0;
            FV(2*inode+1)=-3000;
        }
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
MatrixXd CondenseStiffnessMatrix(vector <Node> NV,MatrixXd KG,MatrixXd FDOF)
{
    int nfreedof=FDOF.sum();
    MatrixXd FKG=MatrixXd::Zero(nfreedof,nfreedof);
    int icounter=-1,jcounter=-1;

    for (int idof=0;idof<2*NV.size();idof++)
    {
        if (FDOF(idof)==1)
        {
            icounter++;
            for (int jdof=0;jdof <2*NV.size();jdof++)
            {
                if (FDOF(jdof)==1)
                {
                    jcounter++;
                    FKG(icounter,jcounter)=KG(idof,jdof);
                }
            }
            jcounter=-1;
        }
    }
    return FKG;
}
MatrixXd Solve_CSE(vector <Node> NV,vector <CSE_Element> EV,MatrixXd FG,MatrixXd KG,MatrixXd FDOF)
{
    MatrixXd UG=MatrixXd::Zero(FDOF.sum(),1);
    MatrixXd UT=MatrixXd::Zero(2*NV.size(),1);
    MatrixXd U=MatrixXd::Zero(6,1);
    UG=KG.inverse()*FG;
    vector <Node> obj;
    int id=-1,dof;

    for (int idof=0;idof<2*NV.size();idof++)
    {
        if (FDOF(idof)==1)
        {
            id++;
            UT(idof)=UG(id);
        }
    }

    for (int ielem=0;ielem<EV.size();ielem++)
    {
        obj=EV[ielem].GetNodeObj();

        for (int inode=0;inode<3;inode++)
        {
            U(2*inode)=UT(2*(obj[inode].GetId()-1));
            U(2*inode+1)=UT(2*(obj[inode].GetId()-1)+1);
        }
        EV[ielem].SetU(U);
        EV[ielem].CalculateStrainTensor();
        EV[ielem].CalculateStressTensor();
        EV[ielem].CalculatePrinStress("plane stress");
        EV[ielem].CalculatePrinStrain("plane stress");
    }

    return UT;
}


MatrixXd Solve_Q4(vector <Node> NV,vector <Q4_Element> EV,MatrixXd FG,MatrixXd KG,MatrixXd FDOF)
{
    MatrixXd UG=MatrixXd::Zero(FDOF.sum(),1);
    MatrixXd UT=MatrixXd::Zero(2*NV.size(),1);
    MatrixXd U=MatrixXd::Zero(8,1);
    UG=KG.inverse()*FG;
    vector <Node> obj;
    int id=-1,dof;

    for (int idof=0;idof<2*NV.size();idof++)
    {
        if (FDOF(idof)==1)
        {
            id++;
            UT(idof)=UG(id);
        }
    }

    for (int ielem=0;ielem<EV.size();ielem++)
    {
        obj=EV[ielem].GetNodalObj();

        for (int inode=0;inode<4;inode++)
        {
            U(2*inode)=UT(2*(obj[inode].GetId()-1));
            U(2*inode+1)=UT(2*(obj[inode].GetId()-1)+1);
        }
        EV[ielem].SetU(U);
        EV[ielem].Get_Ep();
        EV[ielem].Get_Sigma();
        EV[ielem].Set_PSigma("pstress");
        EV[ielem].Set_PStrain("pstress");
    }

    return UT;
}
MatrixXd Solve_Q9(vector <Node> NV,vector <Q4_Element> EV,MatrixXd FG,MatrixXd KG,MatrixXd FDOF)
{
    MatrixXd UG=MatrixXd::Zero(FDOF.sum(),1);
    MatrixXd UT=MatrixXd::Zero(2*NV.size(),1);
    MatrixXd U=MatrixXd::Zero(8,1);
    UG=KG.inverse()*FG;
    vector <Node> obj;
    int id=-1,dof;

    for (int idof=0;idof<2*NV.size();idof++)
    {
        if (FDOF(idof)==1)
        {
            id++;
            UT(idof)=UG(id);
        }
    }

    for (int ielem=0;ielem<EV.size();ielem++)
    {
        obj=EV[ielem].GetNodalObj();

        for (int inode=0;inode<4;inode++)
        {
            U(2*inode)=UT(2*(obj[inode].GetId()-1));
            U(2*inode+1)=UT(2*(obj[inode].GetId()-1)+1);
        }
        EV[ielem].SetU(U);
        EV[ielem].Get_Ep();
        EV[ielem].Get_Sigma();
        EV[ielem].Set_PSigma("pstress");
        EV[ielem].Set_PStrain("pstress");
    }

    return UT;
}
