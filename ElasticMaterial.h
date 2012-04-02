#ifndef ELASTICMATERIAL_H
#define ELASTICMATERIAL_H


class ElasticMaterial
{
    public:
        ElasticMaterial();
        ElasticMaterial(double=2e9,double=0.3,double=8050); //E,v,density in kg/m^3
        void SetE(double);
        double GetE(void);
        void Setv(double);
        double Getv(void);
        void SetGama(double);
        double GetGama(void);
        virtual ~ElasticMaterial();

    protected:

    private:
        double E,v,gama;
};

#endif // ELASTICMATERIAL_H
