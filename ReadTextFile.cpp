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
    Node obj;

    while (!fid.eof())
    {
        fid >> id >> xcord >> ycord >>zcord;
        obj.SetId(id);
        obj.SetCord(xcord,ycord);
        obj.SetForce(0.0,0.0);
        obj.SetRes(xr,yr);
        NV.push_back(obj);
    }

    return NV;
}

vector <CSE_Element> ReadCSEElement(string str,ElasticMaterial mat,double thickness,vector <Node> N)
{
    ifstream fid;
    fid.open(str);

    if (!fid) cout << "Can''t open the CSE element  file" <<endl;
    if (fid) cout <<  "Open the CSE element file" <<endl;

    int id;
    int n1,n2,n3;
    vector <CSE_Element> CSEV;
    CSE_Element obj;

    while (!fid.eof())
    {
        fid >> id >> n1 >> n2 >> n3;
        obj.SetId(id);
        obj.Setth(thickness);
        obj.SetNodeObj(N[n1-1],N[n2-1],N[n3-1]);
        obj.SetMat(mat);
        CSEV.push_back(obj);
    }
    return CSEV;
}

vector <Q4_Element>  ReadQ4Element(string str,ElasticMaterial mat,double thickness,vector <Node> NV)
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
        obj.SetNodalObj(NodeObj);
        obj.SetMat(mat);
        Q4V.push_back(obj);
    }
    return Q4V;
}

vector <Q8_Element>  ReadQ8Element(string str,ElasticMaterial mat,double thickness,vector <Node> NV)
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
        obj.SetNodalObj(NodeObj);
        obj.SetMat(mat);
        Q8V.push_back(obj);
    }
    return Q8V;
}

vector <Q9_Element>  ReadQ9Element(string str,ElasticMaterial mat,double thickness,vector <Node> NV)
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
        obj.SetNodalObj(NodeObj);
        obj.SetMat(mat);
        Q9V.push_back(obj);
    }
    return Q9V;
}
