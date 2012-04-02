#include "ElasticMaterial.h"

ElasticMaterial::ElasticMaterial()
{
    //ctor
}
ElasticMaterial(double Ev,double vv,double g)
{
    SetE(Ev);
    Setv(vv);
    SetGama(g);
}
void ElasticMaterial::SetE(double Ev)
{
    E=Ev;
}
double ElasticMaterial::GetE(void)
{
    return E;
}
void ElasticMaterial::Setv(double vv)
{
    v=vv;
}
double ElasticMaterial::Getv(void)
{
    return v;
}
void ElasticMaterial::SetGama(double g)
{
    gama=g;
}
double ElasticMaterial::GetGama(void)
{
    return gama;
}
ElasticMaterial::~ElasticMaterial()
{
    //dtor
}
