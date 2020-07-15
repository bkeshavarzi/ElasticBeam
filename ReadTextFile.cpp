#include "TextFiles.h"

vector <Node> ReadNodeFile(string str)
{
    ifstream fid;
    fid.open(str);

    if (!fid) cout << "Can''t open the node text file" <<endl;
    if (fid)  cout << "Open the node text file" <<endl;

    int id;
    double xcord,ycord,zcord;
    bool xr=false,yr=false;
    vector <Node> NV;

    while (!fid.eof())
    {
        Node obj;
        fid >> id >> xcord >> ycord >>zcord;
        obj.SetId(id);
        obj.SetCord(xcord,ycord);
        obj.SetForce(0.0,0.0);
        obj.SetRes(xr,yr);
        NV.push_back(obj);
    }

    return NV;
}

vector <CSE_Element> ReadCSEElement(string str,ElasticMaterial & mat,double thickness,vector <Node> N)
{
    ifstream fid;
    fid.open(str);

    if (!fid) cout << "Can''t open the CSE element  file" <<endl;
    if (fid) cout <<  "Open the CSE element file" <<endl;

    int id;
    int n1,n2,n3;
    vector <CSE_Element> CSEV;

    while (!fid.eof())
    {
        CSE_Element obj;
        fid >> id >> n1 >> n2 >> n3;
        obj.SetId(id);
        obj.Setth(thickness);
        obj.SetNodeObj(N[n1-1],N[n2-1],N[n3-1]);
        obj.SetMat(mat);
        obj.SetLocalCord();
        obj.SetDMatrix_PlaneStress();
        CSEV.push_back(obj);
    }
    return CSEV;
}

vector <Q4_Element>  ReadQ4Element(string str,ElasticMaterial & mat,double thickness,vector <Node> NV)
{
    ifstream fid;
    fid.open(str);

    if (!fid) cout << "Can''t open the Q4 element  file" <<endl;
    if (fid) cout <<  "Open the Q4 element  file" <<endl;

    int id;
    int n1,n2,n3,n4;
    vector <Q4_Element> Q4V;
    vector <Node> NodeObj;
    NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);
    Q4_Element obj;

    while (!fid.eof())
    {
        fid >> id >> n1 >> n2 >> n3 >> n4;
        obj.SetId(id);
        obj.Setth(thickness);
        NodeObj[0]=NV[n1-1];NodeObj[1]=NV[n2-1];NodeObj[2]=NV[n3-1];NodeObj[3]=NV[n4-1];
        NodeObj=SortElementNodeQ4(NodeObj);
        obj.SetNodalObj(NodeObj);
        obj.SetMat(mat);
        obj.Setlocalcord();
        obj.SetDMatrix("pstress");
        Q4V.push_back(obj);
    }
    return Q4V;
}

vector <Q8_Element>  ReadQ8Element(string str,ElasticMaterial & mat,double thickness,vector <Node> NV)
{
    ifstream fid;
    fid.open(str);

    if (!fid) cout << "Can''t open the Q8 element  file" <<endl;
    if (fid) cout <<  "Open the Q8 element  file" <<endl;

    int id;
    int n1,n2,n3,n4,n5,n6,n7,n8;
    vector <Q8_Element> Q8V;
    vector <Node> NodeObj;
    NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);
    Q8_Element obj;

    while (!fid.eof())
    {
        fid >> id >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7 >> n8;
        obj.SetId(id);
        obj.Setth(thickness);
        NodeObj[0]=NV[n1-1];NodeObj[1]=NV[n2-1];NodeObj[2]=NV[n3-1];NodeObj[3]=NV[n4-1];NodeObj[4]=NV[n5-1];NodeObj[5]=NV[n6-1];NodeObj[6]=NV[n7-1];NodeObj[7]=NV[n8-1];
        NodeObj=SortElementNodeQ8(NodeObj);
        obj.SetNodalObj(NodeObj);
        obj.SetMat(mat);
        obj.SetLocalCord();
        obj.SetDMatrix("pstress");
        Q8V.push_back(obj);
    }
    return Q8V;
}

vector <Q9_Element>  ReadQ9Element(string str,ElasticMaterial & mat,double thickness,vector <Node> NV)
{
    ifstream fid;
    fid.open(str);

    if (!fid) cout << "Can''t open the Q9 element  file" <<endl;
    if (fid) cout <<  "Open the Q9 element  file" <<endl;

    int id;
    int n1,n2,n3,n4,n5,n6,n7,n8,n9;
    vector <Q9_Element> Q9V;
    vector <Node> NodeObj;
    NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);NodeObj.push_back(NV[0]);
    Q9_Element obj;

    while (!fid.eof())
    {
        fid >> id >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7 >> n8 >> n9;
        obj.SetId(id);
        obj.Setth(thickness);
        NodeObj[0]=NV[n1-1];NodeObj[1]=NV[n2-1];NodeObj[2]=NV[n3-1];NodeObj[3]=NV[n4-1];NodeObj[4]=NV[n5-1];NodeObj[5]=NV[n6-1];NodeObj[6]=NV[n7-1];NodeObj[7]=NV[n8-1];NodeObj[8]=NV[n9-1];
        NodeObj=SortElementNodeQ9(NodeObj);
        obj.SetNodalObj(NodeObj);
        obj.SetMat(mat);
        obj.Setlocalcord();
        obj.SetDMatrix("pstress");
        Q9V.push_back(obj);
    }
    return Q9V;
}

vector <Node> SortElementNodeQ4(vector <Node> Nobj)
{
    int n=Nobj.size();
    MatrixXd mat=MatrixXd::Zero(2,4);
    vector <double> cord;
    for (int inode=0;inode<n;inode++)
    {
        cord=Nobj[inode].GetCord();
        mat(0,inode)=cord[0];
        mat(1,inode)=cord[1];
    }

    double x_min=mat.block(0,0,1,4).minCoeff();
    double x_max=mat.block(0,0,1,4).maxCoeff();
    double y_max=mat.block(1,0,1,4).maxCoeff();
    double y_min=mat.block(1,0,1,4).minCoeff();

    vector <Node> OutNode=Nobj;

    for (int inode=0;inode<n;inode++)
    {
        cord=Nobj[inode].GetCord();
        if ((cord[0]==x_min)&&(cord[1]==y_min)) OutNode[0]=Nobj[inode];
        if ((cord[0]==x_max)&&(cord[1]==y_min)) OutNode[1]=Nobj[inode];
        if ((cord[0]==x_max)&&(cord[1]==y_max)) OutNode[2]=Nobj[inode];
        if ((cord[0]==x_min)&&(cord[1]==y_max)) OutNode[3]=Nobj[inode];
    }

    return OutNode;
}

vector <Node> SortElementNodeQ8(vector <Node> Nobj)
{
    int n=Nobj.size();
    MatrixXd mat=MatrixXd::Zero(2,n);
    vector <double> cord;
    for (int inode=0;inode<n;inode++)
    {
        cord=Nobj[inode].GetCord();
        mat(0,inode)=cord[0];
        mat(1,inode)=cord[1];
    }

    double x_min=mat.block(0,0,1,n).minCoeff();
    double x_max=mat.block(0,0,1,n).maxCoeff();
    double y_max=mat.block(1,0,1,n).maxCoeff();
    double y_min=mat.block(1,0,1,n).minCoeff();
    double x_ave=0.5*(x_min+x_max);
    double y_ave=0.5*(y_min+y_max);

    vector <Node> OutNode=Nobj;

    for (int inode=0;inode<n;inode++)
    {
        cord=Nobj[inode].GetCord();
        if ((cord[0]==x_min)&&(cord[1]==y_min)) OutNode[0]=Nobj[inode];
        if ((cord[0]==x_ave)&&(cord[1]==y_min)) OutNode[1]=Nobj[inode];
        if ((cord[0]==x_max)&&(cord[1]==y_min)) OutNode[2]=Nobj[inode];
        if ((cord[0]==x_max)&&(cord[1]==y_ave)) OutNode[3]=Nobj[inode];
        if ((cord[0]==x_min)&&(cord[1]==y_ave)) OutNode[4]=Nobj[inode];
        if ((cord[0]==x_min)&&(cord[1]==y_max)) OutNode[5]=Nobj[inode];
        if ((cord[0]==x_ave)&&(cord[1]==y_max)) OutNode[6]=Nobj[inode];
        if ((cord[0]==x_max)&&(cord[1]==y_max)) OutNode[7]=Nobj[inode];
    }

    return OutNode;
}

vector <Node> SortElementNodeQ9(vector <Node> Nobj)
{
    int n=Nobj.size();
    MatrixXd mat=MatrixXd::Zero(2,n);
    vector <double> cord;
    for (int inode=0;inode<n;inode++)
    {
        cord=Nobj[inode].GetCord();
        mat(0,inode)=cord[0];
        mat(1,inode)=cord[1];
    }

    double x_min=mat.block(0,0,1,n).minCoeff();
    double x_max=mat.block(0,0,1,n).maxCoeff();
    double y_max=mat.block(1,0,1,n).maxCoeff();
    double y_min=mat.block(1,0,1,n).minCoeff();
    double x_ave=0.5*(x_min+x_max);
    double y_ave=0.5*(y_min+y_max);

    vector <Node> OutNode=Nobj;

    for (int inode=0;inode<n;inode++)
    {
        cord=Nobj[inode].GetCord();
        if ((cord[0]==x_min)&&(cord[1]==y_min)) OutNode[0]=Nobj[inode];
        if ((cord[0]==x_ave)&&(cord[1]==y_min)) OutNode[1]=Nobj[inode];
        if ((cord[0]==x_max)&&(cord[1]==y_min)) OutNode[2]=Nobj[inode];
        if ((cord[0]==x_max)&&(cord[1]==y_ave)) OutNode[3]=Nobj[inode];
        if ((cord[0]==x_ave)&&(cord[1]==y_ave)) OutNode[4]=Nobj[inode];
        if ((cord[0]==x_min)&&(cord[1]==y_ave)) OutNode[5]=Nobj[inode];
        if ((cord[0]==x_min)&&(cord[1]==y_max)) OutNode[6]=Nobj[inode];
        if ((cord[0]==x_ave)&&(cord[1]==y_max)) OutNode[7]=Nobj[inode];
        if ((cord[0]==x_max)&&(cord[1]==y_max)) OutNode[8]=Nobj[inode];
    }

    return OutNode;
}

void WriteOutPutFile_CSE(MatrixXd FU, vector <Node> NV,vector <CSE_Element> EV)
{


    ofstream fidd("Displacement.txt");
    ofstream fids("Stress_Strain.txt");

    vector <double> temp;



    MatrixXd Ep=MatrixXd::Zero(3,1);
    MatrixXd Sigma=MatrixXd::Zero(3,1);
    MatrixXd PSigma=MatrixXd::Zero(3,1);
    MatrixXd PStrain=MatrixXd::Zero(3,1);



    for (int inode=0;inode<NV.size();inode++)
    {
        temp=NV[inode].GetCord();
        fidd << inode << "\t" << temp[0] << "\t" << temp[1] << "\t" << FU(2*inode) << "\t" << FU(2*inode+1) <<endl;
    }

    fidd.close();

    for (int ielem=0;ielem<EV.size();ielem++)
    {
        Ep=EV[ielem].CalculateStrainTensor();
        Sigma=EV[ielem].CalculateStressTensor();
        PSigma=EV[ielem].CalculatePrinStress("plane stress");
        PStrain=EV[ielem].CalculatePrinStrain("plane stress");
        fids <<Ep(0)<<"\t"<<Ep(1)<<"\t"<<Ep(2)<<"\t"<<Sigma(0)<<"\t"<<Sigma(1)<<"\t"<<Sigma(2)<<"\t"<<PSigma(0)<<"\t"<<PSigma(1)<<"\t"<<PSigma(2)<<"\t"<<PStrain(0)<<"\t"<<PStrain(1)<<"\t"<<PStrain(2)<<endl;
    }

    fids.close();
}
