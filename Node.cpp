#include "Node.h"

Node::Node()
{
    //ctor
}
Node::Node(int i,double x,double y,double xf,double yf,bool xr,bool yr)
{
    SetId(i);
    SetCord(x,y);
    SetForce(xf,yf);
    SetRes(xr,yr);
}
void Node::SetId(int i)
{
    id=i;
}
int Node::GetId()
{
    return id;
}
void Node::SetCord(double x,double y)
{
    xc=x;
    yc=y;
}
vector <double> Node:: GetCord(void)
{
    vector <double> temp;
    temp.push_back(xc);
    temp.push_back(yc);
    return temp;
}
void Node::SetForce(double xf,double yf)
{
    fx=xf;
    fy=yf;
}
vector <double> Node:: GetForce(void)
{
    vector <double> temp;
    temp.push_back(fx);
    temp.push_back(fy);
    return temp;
}
void Node::SetRes(bool xr,bool yr)
{
    rx=xr;
    ry=yr;
}
vector <bool>Node:: GetRes(void)
{
    vector <bool> temp;
    temp.push_back(rx);
    temp.push_back(ry);
    return temp;
}

Node::~Node()
{
    //dtor
}
